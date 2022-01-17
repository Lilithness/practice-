
"""
Get descriptions about the variants found between two DNA sequences.
"""

import io
import os
import sys
import json
from typing import Dict, List, Any, Sequence
import argparse
from datetime import datetime
import requests
import pandas as pd
from corvus.cmd import get_cmd_output  # type: ignore

timestamp = datetime.now()
## ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ##
def read_vcf(path: str) -> pd.DataFrame:
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
def symbol_database(symbol: str) -> None:
    """
    Find the species and database for a symbol in a linked external database
    :param: Gene symbol
    :return: None
    """
    url = "https://rest.Ensembl.org/lookup/symbol/homo_sapiens/"  # to get the needed info about the gene chosen to get a sequence from Ensembl
    req = requests.get(url + symbol, headers={"Content-Type": "application/json"})
    req.raise_for_status()
    return req.json()
def get_reference(symbol: str) -> None:
    """
    Get the reference sequence as a FASTA file from Ensembl
    :param: Gene symbol
    :return: None
    """
    interval = f'{symbol_database(symbol)["seq_region_name"]}:{symbol_database(symbol)["start"]}-{symbol_database(symbol)["end"]}:{symbol_database(symbol)["strand"]}'
    url = "https://rest.Ensembl.org/sequence/region/human/"
    req = requests.get(url + interval, headers={"Content-Type": "text/x-fasta"})
    req.raise_for_status()
    return req.text




## ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ##
def main() -> None:
    """
    Check if a DNA sequence has a mutation in the hemoglobin-Beta gene
    :param query: Query file path
            symbol: Gene symbol
            enst: Ensembl transcript id
    :return: None
    """

    parser = argparse.ArgumentParser(description='Get descriptions about the variants found between two DNA sequences')
    parser.add_argument('-q', '--query', type=str, help='Path to the query FASTA file',
                        required=True)  # TODO: implement support for FASTQ
    parser.add_argument('-s', '--symbol', type=str, help='Symbol of the gene of interest', required=True)
    parser.add_argument('-t', '--enst', type=str, help='Ensemble transcript id of interest', required=True)
    args = vars(parser.parse_args())

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width',100) #TODO: change this number to smth more reasonable


    reference_file = f'{args["symbol"]}_reference.fsa'
    with open(reference_file, "w", encoding='ascii') as file:
        textfile = get_reference(args["symbol"])
        file.write(textfile)

    alignment_result = f'{args["symbol"]}{timestamp}_.vcf'  # the name of the alignment file depending on the gene symbol
    output = get_cmd_output(f"bio align {args['query']} {reference_file} --vcf > {alignment_result}")  # check if the colomn option was necessary

    if output["rc"]:
        print("[ stderr ]")
        print(output["stderr"])
        sys.exit(1)

    print("[ stdout ]")
    for line in output["stdout"]:
        print(line)

    print(f"Working with input file: '{args['query']}'")


    decoded = symbol_database(args['symbol'])
    vcf_df = read_vcf(alignment_result)
    vcf_df["POS"] += int(decoded["start"])  # getting the genomic indexing. had to duplicate the same code from fetch seq

    variants: List = []  # list of all possible variants of the sequences in the specific right format
    for row in vcf_df.itertuples():
        region = row.CHROM
        chrom = region.split(":")[2]
        strand = region.split(":")[-1]
        end_pos = row.POS + len(row.REF) - 1
        line = f"{chrom} {row.POS} {end_pos} {row.REF}/{row.ALT} {strand}"  # TODO: check if smth really requires the forward slash downstream -> Crucial for the format of the variant
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
            example["dict_of_transcript_consequences"] = [
                item for item in example['transcript_consequences'] if isinstance(item, dict)
            ]       # change the value of the dictionary to be a dictionary

            ## Extract consequences into a list of dictionaries to be appended to the main variant list datastructure
            example['significant_transcript_consequences'] = [
                item for item in example['dict_of_transcript_consequences'] if item['transcript_id'] == args['enst']  # ENST00000335295
            ][0]       # item is a dict # get only the dictionaries that contain a specific value for transcript id key
            consequences_list.append(example['significant_transcript_consequences'])

            if 'regulatory_feature_consequences' in example:
                example['rf_cons'] = example['regulatory_feature_consequences'][0]  # to correct the format of the output.


            example['Ref'] = example['allele_string'].split("/")[0]
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
        filename = f'{args["symbol"]}_{timestamp}.annotated.csv'
        output_df.to_csv(filename, index=False)
        print(output_df)


## ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ##
if __name__ == "__main__":
    main()
