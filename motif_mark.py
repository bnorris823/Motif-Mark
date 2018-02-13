import argparse
import



def get_parts(line):
    '''Split fasta sequence into pre-intron, exon, and post-intron parts. returned as a list.'''
    pre_intron = ""
    post_intron = ""
    exon = ""
    ex_found = False

    for c in line:
        if c.isupper():
            exon = exon + c
            ex_found = True
        elif ex_found:
            post_intron = post_intron + c
        else:
            pre_intron = pre_intron + c

    return [pre_intron, exon, post_intron]





def fasta_to_dict(file):
    '''Takes a fasta file and makes each entry a key value pair in a dictionary'''
    with open(file)as fh:
        fasta_dict = {}
        header = ""
        seq = ""
        first = 0

        for line in fh:
            line = line.strip()
            if line[0] == ">":
                if not first:
                    fasta_dict[header] = seq
            else:
                seq += line
                first = False

    return fasta_dict




def parse_motifs(motif):
    '''Takes a motif and adds it to a dict with a list of the possible chars at each position as its value.'''
    motif_list = []
    for c in motif:
        motif_list.append(base_dict[c])

    if motif not in motif_dict:
        motif_dict[motif] = motif_list



base_dict = {"A": "Aa", "T":"Tt", "C":"Cc", "G":"Gg", "U":"Uu",
              "R": "AaGg", "Y":"TtCcUu", "S":"CcGg", "W":"AaTtUu",
              "K":"GgTtUu", "M":"AaCc", "B":"CcGgTtUu", "D":"AaGgTtUu",
              "H":"AaCcTtUu", "V":"AaCcGg", "N":"AaTtCcGgUu"}

motif_dict = {}

