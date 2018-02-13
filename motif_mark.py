
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
        count = 0
        for line in fh:
            line = line.strip()
            if line[0] == ">":
                if count > 0:
                    fasta_dict[header] = seq
            else:
                seq = line
                count += 1
    return fasta_dict
