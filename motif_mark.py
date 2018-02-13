



def get_parts(line):
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

def parse_motif(motif):
    motif_list = []
    for c in motif:
        motif_list.append(base_dict[c])

    return motif_list


base_dict = {"A": "Aa", "T":"Tt", "C":"Cc", "G":"Gg", "U":"Uu",
              "R": "AaGg", "Y":"TtCcUu", "S":"CcGg", "W":"AaTtUu",
              "K":"GgTtUu", "M":"AaCc", "B":"CcGgTtUu", "D":"AaGgTtUu",
              "H":"AaCcTtUu", "V":"AaCcGg", "N":"AaTtCcGgUu"}

