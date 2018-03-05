import argparse
import cairo


def get_arguments():
    parser = argparse.ArgumentParser(description="Python script to display binding site motifs around exons")
    parser.add_argument("-f", help="set fasta filename", required=True, type=str)
    parser.add_argument("-m", help="set motif file", required=False, type=bool)

    return parser.parse_args()


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



def get_exon(line):
    '''finds exon in fasta sequence and return it'''
    exon = ""

    for c in line:
        if c.isupper():
            exon = exon + c

    return exon




def fasta_to_dict(file):
    '''Takes a fasta file and makes each entry a key value pair in a dictionary'''
    with open(file)as fh:
        fasta_dict = {}
        header = ""
        seq = ""
        first = True

        for line in fh:
            line = line.strip()

            if line[0] == ">":
                if not first:
                    fasta_dict[header] = get_parts(seq)
                header = line
                seq = ""
            else:
                seq += line
                first = False

    fasta_dict[header] = get_parts(seq)
    return fasta_dict




def parse_motifs(motif_file):
    '''Takes a list of motifs and creates a list of possible chars for each motif.'''
    motif_list = []
    motif_dict = {}
    with open(motif_file)as fh:
        for line in fh:
            line = line.strip()
            line = line.upper()
            motif_list.append(line)

    for motif in motif_list:
        chars = []
        for c in motif:
            chars.append(base_dict[c])
        motif_dict[motif] = chars


    return motif_dict






def get_positions(header, seq):
    positions = {}
    for motif in motif_dict.keys():
        chars = motif_dict[motif]
        window = len(motif)
        for i in range(0,len(seq) - window):
            pattern = seq[i:i + window]

            count = 0
            match = True
            for c in pattern:

                if c not in chars[count]:
                    match = False
                count += 1
            if match:

                if header + "_" + motif in positions:
                    positions[header + "_" + motif].append(i)
                else:
                    positions[header + "_" + motif] = [i]

    return positions



base_dict = {"A": "Aa", "T":"Tt", "C":"Cc", "G":"Gg", "U":"Uu",
              "R": "AaGg", "Y":"TtCcUu", "S":"CcGg", "W":"AaTtUu",
              "K":"GgTtUu", "M":"AaCc", "B":"CcGgTtUu", "D":"AaGgTtUu",
              "H":"AaCcTtUu", "V":"AaCcGg", "N":"AaTtCcGgUu"}


#args = get_arguments()

fasta_file= "test.fa" #args.f
motif_file = "motifs_test.txt" #args.m

def find_all_motifs(fasta_file, motifs_file):
    fasta_dict = fasta_to_dict("test.fa")
    motif_dict = parse_motifs("motifs_test.txt")
    fasta_pos_list = []

    for item in fasta_dict.items():
        seq = "".join(fasta_dict[item[0]])
        pos = get_positions(item[0], seq)
        pos["exon"] = len(get_exon(seq))
        fasta_pos_list.append(pos)

    return fasta_pos_list