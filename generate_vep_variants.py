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


## ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ##
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


## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ##
# TODO: rename this func to reflect the fact that we are getting gene metadata, not sequences
def fetch_seq(symbol: str) -> None:
    """
    Get the DNA sequence of that symbol from Ensembl in a file format
    :param: Gene symbol
    :return: None
    """
    url = "https://rest.Ensembl.org/lookup/symbol/homo_sapiens/"  # to get the needed info about the gene chosen to get a sequence from Ensembl
    req = requests.get(url + symbol, headers={"Content-Type": "application/json"})
    req.raise_for_status()

    decoded = req.json()
    # TODO: abstract the below code away
    # return req.json()

    interval = f'{decoded["seq_region_name"]}:{decoded["start"]}-{decoded["end"]}:{decoded["strand"]}'
    url = "https://rest.Ensembl.org/sequence/region/human/"
    req = requests.get(url + interval, headers={"Content-Type": "text/x-fasta"})
    req.raise_for_status()

    ## TODO: return the processed output and save data to disk outside of the func
    with open("generated_query.fsa", "w", encoding='ascii') as file:
        file.write(req.text)


## ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ##
def main(query, symbol, transcript_id) -> None:
    """
    Check if a DNA sequence has a mutation in the hemoglobin-Beta gene
    :param query: query file path
            symbol: gene symbol
            id: transcript id
    :return: None
    """

    parser = argparse.ArgumentParser(description='Get descriptions about the variants found between two DNA sequences')
    parser.add_argument('-q', '--query', type=str, help='Path to the query FASTA file')  # TODO: implement support for FASTQ
    parser.add_argument('-s', '--symbol', type=str, help='Symbol of the gene of interest')
    parser.add_argument('-t', '--enst', type=str, help='Ensemble transcript id of interest')
    args = parser.parse_args()
    # args = vars(parser.parse_args())

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 1000)  ## TODO: change this number to smth more reasonable

#    query = args[0]  # Expected to be the name of a FASTA file
    # TODO: change 'sample' to 'reference' (silly Alex)
    if not os.path.exists(query):
        fetch_seq(symbol)  # Expected to be the gene symbol
        query = "generated_query.fsa"
    sample_path = "test/Sample.fsa"  # should i make it an input?

    # TODO: interpolate a timestamp into the file name
    alignment_result = f'{symbol}.vcf'  # the name of the alignment file depending on the gene symbol
    output = get_cmd_output(f"bio align {query} {sample_path} --vcf > {alignment_result}")  # check if the colomn option was necessary

    if output["rc"]:
        print("[ stderr ]")
        for line in output["stderr"]:
            print(line)
        sys.exit(1)

    print("[ stdout ]")
    for line in output["stdout"]:
        print(line)

    print(f"Working with input file: '{query}'")

    # TODO: Fetch the region and the genomic coordinate offset ("start") from Ensembl
    url = "https://rest.Ensembl.org/lookup/symbol/homo_sapiens/"  # to get the needed info about the gene chosen to get a sequence from Ensembl
    req = requests.get(url + symbol, headers={"Content-Type": "application/json"})
    req.raise_for_status()
    decoded = req.json()
    
    vcf_df = read_vcf(alignment_result)
    vcf_df["POS"] += int(decoded["start"])  # getting the genomic indexing. had to duplicate the same code from fetch seq

    variants: List = []  # list of all possible variants of the sequences in the specific right format
    for row in vcf_df.itertuples():
        region = row.CHROM
        chrom = region.split(":")[2]
        strand = region.split(":")[-1]
        end_pos = row.POS + len(row.REF) - 1
        line = f"{chrom} {row.POS} {end_pos} {row.REF}/{row.ALT} {strand}"  # TODO: check if smth really requires the forward slash downstream
        variants.append(line)

    jdata = json.dumps({
        "variants": variants
    })  # the json format needed for the variants

    url = "https://rest.Ensembl.org/vep/homo_sapiens/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    if not variants:
        print("The query exactly matches the reference assembly")
    else:
        req = requests.post(url, headers=headers, data=jdata)
        if not req.ok:
            req.raise_for_status()
            sys.exit()

        data: Sequence[Any] = req.json()  # data is a list of many dictionaries which contain the data from Ensembl
        consequences_list: List = []  # we need it to have direct access to the whole value(dict)
        example_dict_list: List = []

        for i, item in enumerate(data):
            example: Dict = data[i]
            example["transcript_consequences"] = [
                item for item in example['transcript_consequences'] if isinstance(item, dict)
            ]       # change the value of the dictionary to be a dictionary
            # TODO: try to avoid re-assigning variables, esp. when coercing to another type

            ## Extract consequences into a list of dictionaries to be appended to the main variant list datastructure
            example['significant_transcript_consequences'] = [
                item for item in example['transcript_consequences'] if item['transcript_id'] == transcript_id  # ENST00000335295
            ][0]       # item is a dict # get only the dictionaries that contain a specific value for transcript id key
            consequences_list.append(example['significant_transcript_consequences'])

            if 'regulatory_feature_consequences' in example:
                example['rf_cons'] = example['regulatory_feature_consequences'][0]  # to correct the format of the output.

            # example['allele_string'] = [char for char in example['allele_string'] if char != '/']
            # example['Ref'] = [example['allele_string'][0]][0]  # add a key to the dict (Reference Nuc)
            example['allele_string'] = example['allele_string'].split("/")[0]
            assembled_dict = {**example, **consequences_list[i]}  # adding data to the dict
            example_dict_list.append(assembled_dict)

        output_df = pd.DataFrame(
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

        output_df.rename(
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

        ## TODO: use query filename instead and add .annotated.csv postfix+extension
        filename = f'{symbol}.annotated.csv'
        output_df.to_csv(filename, index=False)
        print(output_df)


## ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ##
if __name__ == "__main__":
    main()
