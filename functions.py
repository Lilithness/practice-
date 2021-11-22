from typing import List
import info
import requests


def fetch_seq() -> str:
    url = "https://rest.ensembl.org/sequence/region/human/11: 5,225,464-5,229,395:-1"
    r = requests.get(url, headers={"Content-Type": "text/x-fasta"})
    r.raise_for_status()
    with open("SickleCell.txt", "w+") as file:
        file.write(r.text)
    with open("SickleCell.txt", "r") as file:
        dna = ''.join([i for i in file][1:])  # starting from index 1 to skip the file's header
    return dna


def transcription(seq: str) -> str:
    """
    Replace all "T" characters with a single "U".
    :param seq: a string of letters representing DNA nucleotides
    :return: a string with all "T"s replaced by single "U"
    """
    return seq.replace("T", "U")


def complementary_DNA(seq: str) -> str:
    """
    Make a reverse complement of the given nucleotide sequence.
    :param seq: a DNA sequence in the IUPAC nucleotide code representation
    :return: reverse-complementary DNA sequence
    """
    return ''.join([info.dna_complement_strand[nuc] for nuc in seq])


def translate(seq: str, init_pos=0) -> List[str]:
    """
    Map codons to amino-acids, using IUPAC representations.

    :param seq: a string of letters representing DNA nucleotides
    :param init_pos: index of starting nucleotide
    :return: a string representing the amino-acid sequence
    """
    return [
        info.table[seq[pos:pos+3]] for pos in range(init_pos, len(seq)-2, 3)
            ]


def gen_reading_frames(seq: str) -> List[List[str]]:
    """
    Given a nucleotide sequence, return a list of protein candidates,
    ...corresponding to all 6 possible reading frames.

    :param seq: a DNA sequence in the IUPAC nucleotide code representation
    :return: a list of protein products in IUPAC amino acid representation
    """
    return [
        translate(seq, 0),
        translate(seq, 1),
        translate(seq, 2),
        translate(complementary_DNA(seq), 0),
        translate(complementary_DNA(seq), 1),
        translate(complementary_DNA(seq), 2)
    ]


def proteins_from_rf(aa_seq) -> List[str]:
    """
    Compute all possible proteins. return a list of possible proteins
    :param aa_seq: a string of aminoacids
    :return: a list of all possible proteins
    """
    current_protein: List[str] = []
    proteins: List[str] = []
    for aa in aa_seq:
        for a in aa:
            if a == "_":  # indicates stop codon
                if current_protein:
                    for p in current_protein:
                        proteins.append(p)
                    current_protein = []
            else:
                if a == "M":  # indicates start codon
                    current_protein.append("")
                for i in range(len(current_protein)):
                    current_protein[i] += a
    return proteins


def all_proteins(seq: str, start=0, end=0, ordered=False) -> List[str]:
    """
    Compute possible proteins for 3 open reading frames
    :param: a string representing the DNA sequence
    :return: list of all possible proteins
    """
    if end > start:
        rfs = gen_reading_frames(seq[start:end])
    else:
        rfs = gen_reading_frames(seq)
    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
    if ordered:
        return sorted(res, key=len, reverse=True)
    return res
