from info import *

def transcription(seq: str) -> str:
    """
    Replace all "T" characters with a single "U".
    :param seq: a string of letters representing DNA nucleotides
    :return: a string with all "T"s replaced by single "U"
    """
    return seq.replace("T", "U")

def reverse_compliment(seq):
    return ''.join([dna_reverse_compliment[nuc] for nuc in seq])

def translate(seq, init_pos=0):
    return [table[seq[pos:pos+3]] for pos in range (init_pos, len(seq)-2,3)]

def gen_reading_frames(seq): #6 reading frames
    frames = []
    frames.append(translate(seq, 0))
    frames.append(translate(seq, 1))
    frames.append(translate(seq, 2))
    frames.append(translate(reverse_compliment(seq), 0))
    frames.append(translate(reverse_compliment(seq), 1))
    frames.append(translate(reverse_compliment(seq), 2))
    return frames

def proteins_from_rf(aa_seq):
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

def all_proteins_from_orfs (seq, start=0,end=0,ordered = False):
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
