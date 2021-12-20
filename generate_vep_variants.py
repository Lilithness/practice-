import pandas as pd
from corvus.cmd import get_cmd_output
import vcf
import io
import os
import requests
import sys
import json
import pprint

def read_vcf(path):
    """
    Read a VCF file of alligned sequences
    :param path: VCF file path
    :return: Data frame of all data    
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
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
    r = requests.get(url, headers={"Content-Type": "text/x-fasta"})
    r.raise_for_status()
    with open("SickleCell.txt", "w+") as file:
        file.write(r.text)

def main(args: List[str]) -> None:
    """
    Check if a DNA sequence has a mutation in the hemoglobin-Beta gene
    :param args: list that contains the name of the file
    :return: None
    """

    query = args[0]
    if os.path.exists(query):
        try:
          sample_path = "/home/lilith/input/sample.fsa"
          get_cmd_output(f"bio align {query} {sample_path2} --vcf | column >alignment_result.vcf")

        except FileNotFoundError:
            print(f"Sorry, file not found: '{query}'")
    else:
        dna = fetch_seq()

    print(f"Working with input file: '{query}'")

    pd.set_option('display.max_columns', 60)
    pd.set_option('display.max_rows',None)
    vcf_df = read_vcf("alignment_result.vcf")
    vcf_df["POS"] = vcf_df["POS"] +5225464
    variants = []
    for row in vcf_df.itertuples():
        region = row.CHROM
        chrom = region.split(":")[2]
        strand = region.split(":")[-1]
        ref = row.REF
        pos = row.POS
        end_pos= pos + len(ref)-1
        alt = row.ALT
        line = f"{chrom} {pos} {end_pos} {ref}/{alt} {strand}"
        variants.append(line)
        data = json.dumps({
        "variants": variants
        })
        server = "https://rest.ensembl.org"
        ext = "/vep/homo_sapiens/region"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        if data == '{"variants": []}':
            print("The sequences match")
        else:
            r = requests.post(server+ext, headers=headers, data= data)
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            decoded = r.json()
            for i in decoded:
                pretty_print_json = pprint.pformat(i)
                print(pretty_print_json)
if __name__ == "__main__":
    main(sys.argv[1:])
