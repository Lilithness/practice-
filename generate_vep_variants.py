"""
Get descriptions about the variants found between two DNA sequences.
"""

import io
import os
import sys
import json
from typing import Dict, List, Any, Sequence
import argparse
import requests
import pandas as pd
from corvus.cmd import get_cmd_output  # type: ignore


parser = argparse.ArgumentParser(description='Get descriptions about the variants found between two DNA sequences')
parser.add_argument('query', type=str, help='The path of the query file')
parser.add_argument('symbol', type=str, help='The symbol of the gene of interest')
parser.add_argument('transcript_id', type=str, help='The transcript id of interest ')
args = parser.parse_args()

# TODO: annotate the return value


def read_vcf(path: str):
    """
    Read a VCF file
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


def fetch_seq(symbol: str) -> None:
    """
    Get the DNA sequence of that symbol from Ensembl in a file formate
    :param: Gene symbol
    :return: None
    """
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/"  # to get the needed info about the gene chosen to get a sequence from ensembl
    req = requests.get(url + symbol, headers={"Content-Type": "application/json"})
    req.raise_for_status()
    decoded = req.json()
    interval = f'{decoded["seq_region_name"]}:{decoded["start"]}-{decoded["end"]}:{decoded["strand"]}'
    url = "https://rest.ensembl.org/sequence/region/human/"
    req = requests.get(url + interval, headers={"Content-Type": "text/x-fasta"})
    req.raise_for_status()
    with open("generated_query.fsa", "w", encoding='ascii') as file:
        file.write(req.text)


def main(query, symbol, transcript_id) -> None:
    """
    Check if a DNA sequence has a mutation in the hemoglobin-Beta gene
    :param query: query file path
            symbol: gene symbol
            id: transcript id
    :return: None
    """

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 1000)

#    query = args[0]  # Expected to be the name of a FASTA file
    if not os.path.exists(query):
        fetch_seq(symbol)  # Expected to be the gene symbol
        query = "generated_query.fsa"
    sample_path = "test/Sample.fsa"  # should i make it an input?
    alignment_result = f'alignment_result_{symbol}.vcf'  # the name of the allignment file depending on the gene symbol
    get_cmd_output(f"bio align {query} {sample_path} --vcf > {alignment_result}")  # check if the colomn option was necessary

    print(f"Working with input file: '{query}'")

    # TODO: Fetch the region and the genomic coordinate offset ("start") from Ensembl
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/"  # to get the needed info about the gene chosen to get a sequence from ensembl
    req = requests.get(url + symbol, headers={"Content-Type": "application/json"})
    req.raise_for_status()
    decoded = req.json()
    strpos = int(decoded["start"])
    vcf_df = read_vcf(alignment_result)
    vcf_df["POS"] = vcf_df["POS"] + strpos  # getting the genomic indexing. had to duplicate the same code from fetch seq

    variants: List = []  # list of all possible variants of the sequences in the specific right format
    for row in vcf_df.itertuples():
        region = row.CHROM
        chrom = region.split(":")[2]
        strand = region.split(":")[-1]
        end_pos = row.POS + len(row.REF) - 1
        line = f"{chrom} {row.POS} {end_pos} {row.REF}/{row.ALT} {strand}"
        variants.append(line)

    jdata = json.dumps({
        "variants": variants
    })  # the json format needed for the variants

    url = "https://rest.ensembl.org/vep/homo_sapiens/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    if not variants:
        print("The query exactly matches the reference assembly")
    else:
        req = requests.post(url, headers=headers, data=jdata)
        if not req.ok:
            req.raise_for_status()
            sys.exit()

        data: Sequence[Any] = req.json()  # data is a list of many dictionaries which contain the data from ensembl
        consequences_list: List = []  # we need it to have direct access to the whole value(dict)
        assembled_dict: Any = {}
        example_dict_list: List = []
        for i, item in enumerate(data):

            example: Dict = data[i]  # example is a dictionary
            example["transcript_consequences"] = [
                item for item in example['transcript_consequences'] if isinstance(item, dict)
            ]       # change the value of thed dictionary to be a dictionary
            example['significant_transcript_consequences'] = [
                item for item in example['transcript_consequences'] if item['transcript_id'] == transcript_id  # ENST00000335295
            ][0]       # item is a dict # get only the dictionaries that contain a specific value for transcript id key

            consequences_list.append(example['significant_transcript_consequences'])
            # we need it to have direct access to the whole value(dict) to add it then to the original as keys and values
            # put the dictionaries wanted in a list\

            if not example.keys == 'regulatory_feature_consequences':
                pass
            else:
                example['rf_cons'] = example['regulatory_feature_consequences'][0]  # to correct the format of the output.

            example['allele_string'] = [char for char in example['allele_string'] if char != '/']
            example['Ref'] = [example['allele_string'][0]][0]  # add a key to the dict (Reference Nuc)
            assembled_dict = {**example, **consequences_list[i]}  # adding data to the dict
            example_dict_list.append(assembled_dict)

        new = pd.DataFrame(
            example_dict_list,
            columns=[
               'gene_symbol',
               'seq_region_name',
               'start',
               'end',
               'Ref',
               'variant_allele',
               'strand',
               'biotype',
               'impact'
            ])

        new.rename(
            columns={
                'seq_region_name': 'Chromosome',
                'variant_allele': 'Alt',
                'gene_symbol': 'Gene',
                'start': 'Start',
                'end': 'End',
                'strand': 'Strand',
                'biotype': 'Biotype',
                'impact': 'Impact'
            },
            inplace=True
        )

        filename = f'data_file_{symbol}.csv'
        new.to_csv(filename, index=False)
        print(new)


if __name__ == "__main__":
    main(args.query, args.symbol, args.transcript_id)
