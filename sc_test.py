#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
from functions import fetch_seq
from functions import all_proteins
import requests
from requests.exceptions import HTTPError
import sys
from typing import List
import re
import os.path


get_cmd_output("bio align /home/lilith/input/sample.fsa /home/lilith/input/query.fsa --vcf | column >alignment_result.vcf")
#Run a shell command via subprocess and return exit code and stdout/stderr.

pd.set_option('display.max_columns', 60)
pd.set_option('display.max_rows',None)
#Set pandas options

url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/HBB"
try:
    r = requests.get(url, headers={"Content-Type": "application/json"})
    r.raise_for_status()
except HTTPError as http_err:
    print(f'HTTP error occurred: {http_err}')
except Exception as err:
    print(f'Other error occurred: {err}')
else:
    decoded = r.json()
    print(f'interval: {decoded["seq_region_name"]}: {decoded["start"]}-{decoded["end"]}')
    print(f'type:{decoded["biotype"]} {decoded["object_type"]}')
    print(f'Name of the gene:{decoded["display_name"]}')


def main(args: List[str]) -> None:
    """
    Check if a DNA sequence has a mutation in the hemoglobin-Beta gene
    :param args: list that contains the name of the file
    :return: None
    """

    filename = args[0]
    if os.path.exists(filename):
        try:
            with open(filename, "r") as file:
                dna = file.read().strip().upper()
                re.sub(r"[\n \r]", "", dna)
                dna = dna.replace(" ", "")
        except FileNotFoundError:
            print(f"Sorry, file not found: '{filename}'")
    else:
        dna = fetch_seq()
        dna = dna.replace("\n", "")

    print(f"Working with input file: '{filename}'")

    all_proteins_expected = []
    # Make a list of protein sequences
    for prot in all_proteins(dna, 0, 0, True):
        all_proteins_expected.append(prot)
    for i in all_proteins_expected:
        # Check if the given string is in the list

        if i == "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLVSRLQDRFKETNRNWACGDREDSWVSDRH":
            print("the protein's index in the list is number ", all_proteins_expected.index(i))
            print("no sickle cell anemia")

        if i == "MVHLTPVEKSAVTALWGKVNVDEVGGEALGRLVSRLQDRFKETNRNWACGDREDSWVSDRH":
            print("the protein's index in the list is number ", all_proteins_expected.index(i))
            print("sickle cell anemia positive")


if __name__ == "__main__":
    main(sys.argv[1:])
