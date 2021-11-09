from typing import List
import info


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
    return ''.join([info.dna_compliment_strand[nuc] for nuc in seq])
    

def translate(seq: str, init_pos=0) -> List[str]:
    """
    Map nucleotide triplets (codons) to amino-acids, using IUPAC single-character representations.

    :param seq: a string of letters representing DNA nucleotides
    :param init_pos: ---
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

def proteins_from_rf(aa_seq:List[List[str]]) -> List[str]:
    """
    Compute all possible proteins in an aminoacid seq and return a list of possible proteins
    :param aa_seq: a string of aminoacids
    :return: a list of all possible proteins
    """
    current_protein = []
    proteins = []
    for aa in aa_seq:
        for a in aa:
            if a == "_": #indicates stop codon
                if current_protein:
                    for p in current_protein:
                        proteins.append(p)
                    current_protein = []
            else:
                if a =="M": #indicates start codon
                    current_protein.append("")
                for i in range(len(current_protein)):
                    current_protein[i] += a
    return proteins

def all_proteins_from_orfs(seq: str, start=0, end=0, ordered=False) -> List[str]:
    """
    Compute all possible proteins for all open reading frames
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
