#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
from functions import *
import requests
from requests.exceptions import HTTPError
import sys
import re
from typing import List

dna = fetch_Seq()
dna = dna.upper()
dna = dna.replace("\n", "")
dna = dna.replace("\r", "")
dna = dna.replace(" ", "")

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
    print(f'type:{decoded["biotype"]} {decoded["object_type"]},Name of the gene:{decoded["display_name"]}')

def main(args: List[str]) -> None:
    """
    Check if a DNA sequence has a mutation in the hemoglobin-Beta gene
    :param args: list that contains the name of the file
    :return: None
    """
    filename = args[0]

    print(f"Working with input file: '{filename}'")

    try:
        with open(filename, "r") as file:
            dna = file.read().strip().upper()
            # Look up replacing multiple characters
            re.sub(r"[\n \r]", "", dna)
            dna = dna.replace(" ", "")
    except FileNotFoundError:
        print(f"Sorry, file not found: '{filename}'")

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
    