"""
Get descriptions about the variants found between two DNA sequences
"""

import io
import os
import sys
import json
import pprint
import requests
import pandas as pd
from corvus.cmd import get_cmd_output # type: ignore


def read_vcf(path):
    """
    Read a VCF file of alligned sequences
    :param path: VCF file path
    :return: Data frame of all data
    """
    with open(path, 'r',encoding='ascii') as file:
        lines = [word for word in file if not word.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def fetch_seq():
    """
    Gets the DNA sequence of interest form ensembl
    :return: A fasta file
    """
    url = "https://rest.ensembl.org/sequence/region/human/11: 5,225,464-5,229,395:-1"
    req = requests.get(url, headers={"Content-Type": "text/x-fasta"})
    req.raise_for_status()
    with open("SickleCell.txt", "w+", encoding='ascii') as file:
        file.write(req.text)


def main(args) -> None:
    """
    Check if a DNA sequence has a mutation in the hemoglobin-Beta gene
    :param args: list that contains the name of the file
    :return: None
    """

    query = args[0]
    if os.path.exists(query):
        try:
            sample_path = "/home/lilith/input/sample.fsa"
            get_cmd_output(f"bio align {query} {sample_path} --vcf | column >alignment_result.vcf")

        except FileNotFoundError:
            print(f"Sorry, file not found: '{query}'")
    else:
        fetch_seq()
        sample_path = "/home/lilith/input/sample.fsa"
        query = "SickleCell.txt"
        get_cmd_output(f"bio align {query} {sample_path} --vcf | column >alignment_result.vcf")

    print(f"Working with input file: '{query}'")

    pd.set_option('display.max_columns', 60)
    pd.set_option('display.max_rows', None)
    vcf_df = read_vcf("alignment_result.vcf")
    vcf_df["POS"] = vcf_df["POS"] + 5225464
    variants = []
    for row in vcf_df.itertuples():
        region = row.CHROM
        chrom = region.split(":")[2]
        strand = region.split(":")[-1]
        ref = row.REF
        pos = row.POS
        end_pos = pos + len(ref)-1
        alt = row.ALT
        line = f"{chrom} {pos} {end_pos} {ref}/{alt} {strand}"
        variants.append(line)
    data = json.dumps({
        "variants": variants
    })
    url = "https://rest.ensembl.org/vep/homo_sapiens/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    if data == '{"variants": []}':
        print("The sequences match")
    else:
        req = requests.post(url, headers=headers, data=data)
        if not req.ok:
            req.raise_for_status()
            sys.exit()
        decoded = req.json()
        for i in decoded:
            pretty_print_json = pprint.pformat(i)
            print(pretty_print_json)


if __name__ == "__main__":
    main(sys.argv[1:])
