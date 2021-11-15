from info import *
import requests

def fetch_Seq(None):->str
     """
    Fetch a DNA sequence from Ensembl
    :return: a string that represents the DNA
    """
    url = "https://rest.ensembl.org/sequence/region/human/11: 5,225,464-5,229,395:-1"
    r = requests.get(url, headers={"Content-Type": "text/x-fasta"})
    r.raise_for_status()
    with open("SickleCell.txt", "w+") as file:
        sequence = r.text
        file.write(sequence)
    with open("SickleCell.txt", "r") as file:
        temp = []
        for i in file:
            temp.append(i)
    dna = ''.join(temp[1:])
    return dna
    
    
def transcription(seq: str) -> str:
    """
    Replace all "T" characters with a single "U".
    :param seq: a string of letters representing DNA nucleotides
    :return: a string with all "T"s replaced by single "U"
    """
    return seq.replace("T", "U")

def reverse_compliment(seq: str)-> str:
    """
    uses the dna_revers_compliment dictionary to replace every nucleotide in the DNA sequence with the complimentary one, adding them to a list.
    and then joining all the variables in this list togeather to return a string
    :param seq: a string of letters representing DNA nucleotides
    :return: a string representing the complimentary DNA
    """
    return ''.join([dna_reverse_compliment[nuc] for nuc in seq])
    

def translate(seq: str, init_pos=0)-> str:
    """
    this functions uses the dictionary named table to generate an amino acid sequence for every codon (3 nucleotides) in the given DNA sequence
    it goes from the start of the sequence, detecting every 3 nucleotides at a time, and changing hem into Anminoacids depending on the table already given.
    :param seq: a string of letters representing DNA nucleotides
    :return: a string representing the aminoacid sequence
    """
    return [table[seq[pos:pos+3]] for pos in range (init_pos, len(seq)-2,3)]

def gen_reading_frames(seq: str) -> list:
    """
    uses the translation function to turn the DNA sequence into amino acid sequence and then add it to a list.
    this process repeates 3 times on the exact DNA sequence and 3 times on the cDNA to cover all possible reading frames
    :param seq: a string of letters representing DNA nucleotides
    :return: a list that contains 6 possible aminoacid sequences from all reading frames
    """
    frames = []
    frames.append(translate(seq, 0))
    frames.append(translate(seq, 1))
    frames.append(translate(seq, 2))
    frames.append(translate(reverse_compliment(seq), 0))
    frames.append(translate(reverse_compliment(seq), 1))
    frames.append(translate(reverse_compliment(seq), 2))
    return frames

def proteins_from_rf(aa_seq:str)->list:
    """
    it examins every aminoacid sequence given, starts collecting the amino acid sequences in a list called current protein when it encounters an M which represent a start codon.
    this list might contain many aminoacids sequences that the function is still adding to until it encounters a stop codon,
    then every sequence in this lists moves to a list called proteins which will contain all the possible proteins from all the possible DNA reading frames.
    :param aa_seq: a string of aminoacids
    :return: a list of all possible proteins
    """
    current_protein = []
    proteins = []
    for aa in aa_seq:
        if aa == "_": #stop codon
            if current_protein:
                for p in current_protein:
                    proteins.append(p)
                current_protein = []
        else:
            if aa =="M": # start codon
                current_protein.append("")
            for i in range(len(current_protein)):
                current_protein[i]+=aa
    return proteins

def all_proteins_from_orfs (seq: str, start=0,end=0,ordered = False)->list:
    """
    it uses the gen_reading_frames function (so we get aminoacid sequences) on the given sequence and then the proteins_from_rf function to get all proteins,
    adding them to the res list then sorting them -if needed- depending on their length in ascending order.
    :param seq: a string of letters representing the DNA nucleotides
    :return: a list with all proteins
    """
    if end>start:
        rfs = gen_reading_frames(seq[start:end])
    else:
        rfs = gen_reading_frames(seq)
    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
    if ordered:
        return sorted(res,key=len,reverse=True)
    return res
