import sys
from functions import *

from typing import List
def main(args: List[str]) -> None:
    """
   ---
    :param args: list that contains the name of the file
    :return: None
    """
    filename = args[0]

    print(f"Working with input file: '{filename}'")

    try:
        with open(filename, "r") as file:
            dna = file.read().strip()
            # Look up replacing multiple characters
            dna = dna.upper().replace("\n", "").replace("\r", "").replace(" ", "")
    except FileNotFoundError:
        print(f"Sorry, file not found: '{filename}'")
    
    all_proteins_expected = []
    # Make a list of protein sequences
    for prot in all_proteins_from_orfs(dna, 0, 0, True):
        all_proteins_expected.append(prot)

    print(all_proteins_expected)
    for i in all_proteins_expected:
        # Check if the given string is in the list

        if i == "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLVSRLQDRFKETNRNWACGDREDSWVSDRH":
            print("the biological protein's index in the list is number ", all_proteins_expected.index(i),
                    "\nno sickle cell anemia")
    
    
        if i == "MVHLTPVEKSAVTALWGKVNVDEVGGEALGRLVSRLQDRFKETNRNWACGDREDSWVSDRH":
            print(all_proteins_expected.index(i))
            print("sickle cell anemia positive")


if __name__ == "__main__":
    main(sys.argv[1:])
