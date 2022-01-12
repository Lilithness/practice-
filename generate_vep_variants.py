"""
Get descriptions about the variants found between two DNA sequences
"""

import io
import os
import sys
import json
from typing import Dict, List, Any, Sequence
import requests
import pandas as pd
from corvus.cmd import get_cmd_output  # type: ignore


def read_vcf(path):
    """
    Read a VCF file of alligned sequences
    :param path: VCF file path
    :return: Data frame of all data
    """
    with open(path, 'r', encoding='ascii') as file:
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
    jdata = json.dumps({
        "variants": variants
    })
    url = "https://rest.ensembl.org/vep/homo_sapiens/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    if jdata == '{"variants": []}':
        print("The sequences match")
    else:
        req = requests.post(url, headers=headers, data=jdata)
        if not req.ok:
            req.raise_for_status()
            sys.exit()
        decoded = req.json()
        if not os.path.exists("datafile.json"):
            with open("datafile.json", "a", encoding='ascii') as file:
                file.write(json.dumps(decoded, ensure_ascii=False, indent=4))
        with open("datafile.json", 'r', encoding='ascii') as file:
            data: Sequence[dict] = json.load(file)
        consequences_list = []
        for i, item in enumerate(data):
            example: Any = data[i]
            # example['transcript_consequences']: Union[str, slice] = []

            example["transcript_consequences"] = [
                item for item in example['transcript_consequences'] if isinstance(item, dict)
            ]

            example['significant_transcript_consequences'] = [item for item in example['transcript_consequences'] if item['transcript_id'] == 'ENST00000335295'][0]
            consequences_list.append(example['significant_transcript_consequences'])
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        reglist: List[List[Dict[str, Any]]] = []  # list of lists of one dictionary
        ex_dict_list = []
        for i, item in enumerate(data):
            newexample = data[i]
            # list(newexample['regulatory_feature_consequences'])
            reglist.append(newexample['regulatory_feature_consequences'])
            newexample['allele_string'] = [char for char in newexample['allele_string'] if char != '/']
            newexample['Ref'] = [newexample['allele_string'][0]][0]
            assembled_dict = {**newexample, **consequences_list[i], **reglist[i][0]}
            ex_dict_list.append(assembled_dict)
        new = pd.DataFrame(ex_dict_list,
                           columns=['gene_symbol', 'seq_region_name', 'start', 'end', 'Ref', 'variant_allele',  'strand',
                                    'biotype', 'impact'])
        new.rename(
            columns={'seq_region_name': 'Chromosome', 'variant_allele': 'Alt', 'gene_symbol': 'Gene', 'start': 'Start',
                     'end': 'End', 'strand': 'Strand', 'biotype': 'Biotype', 'impact': 'Impact'}, inplace=True)
        new.to_csv('datafile11_new.csv', index=False)
        pd.set_option('display.max_columns', 60)
        print(new)


if __name__ == "__main__":
    main(sys.argv[1:])
