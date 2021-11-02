import sys
import functions

def main(args:list) -> none:
    """
    assign the name of the file
    :param: list that contains the name f the file
    :return: none
    """
    filename = args[0]

with open ("filename", "r") as file:
    """
    open the file that has the DNA sequence and process it
    """
    contents = file.read()
    contents = contents.strip()
    dna = ''.join(contents)
    dna = dna.upper()
    dna = dna.replace("\n", "")
    dna = dna.replace("\r", "")
    dna = dna.replace(" ", "")
    
try:
    """
    try to open the file, send an error massage if it wasn't found
    """
    filename = open(dna, "r")
except FileNotFoundError:
    print("sorry, file not found")

all_protiens_expected = []
"""
make a list of protein sequences
"""
for prot in all_proteins_from_orfs(dna:str, 0, 0, True) -> list[list[str]]:
    all_protiens_expected.append(prot)

for i in all_protiens_expected:
    """
    check if this string given, is in the list
    """
    if i == "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLVSRLQDRFKETNRNWACGDREDSWVSDRH":
        print("the biological protein's index in the list is number ", all_protiens_expected.index(i),
                "\nno sickle cell anemia")


    if i == "MVHLTPVEKSAVTALWGKVNVDEVGGEALGRLVSRLQDRFKETNRNWACGDREDSWVSDRH":
        print(all_protiens_expected.index(i))
        print("sickle cell anemia positive")

if __name__ == "__main__":
    main(sys.argv[1:])
