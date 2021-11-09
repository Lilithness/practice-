#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

from functions import *
import requests
from pprint import pprint
from requests.exceptions import HTTPError

sequence = input("add the file name of the sequence") + ".txt"
sequence = sequence.strip()

try:
    f = open(sequence, "r")
except FileNotFoundError:
    print("sorry, file not found")
    quit()
dna = f.read()
dna = dna.upper()
dna = dna.replace("\n", "")
dna = dna.replace("\r", "")
dna = dna.replace(" ", "")
def main():
    print("What do you want to test the DNA sequence for?")
    print("For Sickle cell anemia enter 1"
          "\nTo exit enter 0")
    return
main()

all_protiens_expected = []
for prot in all_proteins_from_orfs(dna, 0, 0, True):
    all_protiens_expected.append(prot)
#print(all_protiens_expected)
type = int(input())

class Data:
    def __init__(self,disease,mutation_type):
        self.disease = disease
        self.mutation_types = mutation_type
    def full_data(self):
        return '{} {}'.format(self.disease, self.mutation_types)
def check_data(type):
    Q= input("do you want to check your data? Y or N")
    if (Q=="y") or (Q=="Y"):
        if type == 1:
            print(data_1.full_data())

while (type!=0):

    if type == 1:
        for i in all_protiens_expected:
            if i == "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLVSRLQDRFKETNRNWACGDREDSWVSDRH":
                print("the biological protein's index in the list is number ", all_protiens_expected.index(i),
                      "\nno sickle cell anemia")
                data_1 = Data('normal,',' no_mutation')
                check_data(1)

            if i == "MVHLTPVEKSAVTALWGKVNVDEVGGEALGRLVSRLQDRFKETNRNWACGDREDSWVSDRH":
                print(all_protiens_expected.index(i))
                print("sickle cell anemia positive")
                more_info = input("would you like to know more info about this genetic mutation? (Y) or (N)")
                if (more_info == "Y" or more_info == "y"):
                    print(sickle_cell_info)
                data_1 = Data('sickle cell anemia positive,', 'point mutation')
                check_data(1)
    type = 0

pass
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
