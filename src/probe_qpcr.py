fpri = "GGGGAACTTCTCCTGCTAGAAT"
rpri = "CAGACATTTTGCTCTCAAGCTG"
probe = "TTGCTGCTGCTTGACAGATT"

import regex


def read_primers(file_name: str) -> dict:
    """
    Reads primers from file.
    Returns a nested dictionary.
    """
    primer_dict = {}
    with open(file_name) as fin:
        name, seq, direction = None, "", None
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    # print('\tName: {}\n\tDirection: {}\n\tSequence: {}'.format(name, direction, seq))
                    primer_dict[name][direction] = seq
                *name, direction = line.split(".")
                name[0] = name[0][1:]
                name[-1] = name[-1][:-2]
                name = '.'.join(name)
                seq = ""
                # print('\tName: {}\n\tDirection: {}\n\tSequence: {}'.format(name, direction, seq))
                if not primer_dict.get(name):
                    primer_dict[name] = {'F': '', 'R': '', 'P': ''}
            else:
                seq = line
        if name:
            primer_dict[name][direction] = seq

    return primer_dict


def reverse_complement(nuc_sequence: str):
    """
    Returns the reverse complement of a nucleotide sequence.
    >>> reverse_complement('ACGT')
    'ACGT'
    >>> reverse_complement('ATCGTGCTGCTGTCGTCAAGAC')
    'GTCTTGACGACAGCAGCACGAT'
    >>> reverse_complement('TGCTAGCATCGAGTCGATCGATATATTTAGCATCAGCATT')
    'AATGCTGATGCTAAATATATCGATCGACTCGATGCTAGCA'
     """
    complements = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A"
    }
    rev_seq = "".join([complements[s] for s in nuc_sequence[::-1]])
    return rev_seq


def read_sequences(file_name: str):
    """
	Sequence iterable 
	"""
    with open(file_name) as fin:
        name, seq = None, []
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield name, ''.join(seq)
                name, seq = line, []
            else:
                seq.append(line)
        if name:
            yield name, ''.join(seq)


def get_qpcr_hits(primer_name: str, forward_primer: str, reverse_primer: str, probe: str,
                  target_name: str, target_sequence: str):
    """
	Returns the TaqMan product.
	"""
    import re
    # for match in [match for match in re.finditer("({}).*({})".format(
    #         forward_primer, reverse_complement(reverse_primer)), sequence)]:
    #     product = sequence[match.start():match.end()]
    #     for product_match in [match for match in re.finditer(probe, product)]:
    primer_pattern = "({}).*({}).*({})".format( forward_primer, probe, reverse_complement(reverse_primer))
    for match in [match for match in re.finditer(primer_pattern, target_sequence)]:
        product = target_sequence[match.start():match.end()]
        print("Primer set: {}".format(primer_name))
        print("Target sequence: {}".format(target_name))
        print("Successful hit at {}".format(match.start()))
        print("Product: {}".format(product))
        print("\tStart: {}\n\tEnd: {}\n\tLength: {}".format(
            match.start(), match.end(), match.end() - match.start()))

def replace_ambiguity_codes(sequence: str):
    import sys
    codes = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "S": "[GC]", "B": "[CGT]", "Y": "[CT]", "W": "[AT]",
             "D": "[AGT]", "K": "[GT]", "N": "[ACGT]", "H": "[ACT]", "M": "[AC]", "V": "[ACG]", "X": "[ACGT]"}
    regexlist = []
    for base in sequence:
        if not base in codes:
            print("Unidentified nucleotide code in primer:", base)
            sys.exit()
        else:
            regexlist.append(codes[base])

    return ''.join(regexlist)


if __name__ == "__main__":
    primer_dict = read_primers("../data/collected_primers.N.formatted.fa")
    for target_name, target_sequence in read_sequences("../data/test_sequence.fa"):
        for primer_set in primer_dict:
            get_qpcr_hits(primer_set,
                          primer_dict[primer_set]['F'],
                          primer_dict[primer_set]['R'],
                          primer_dict[primer_set]['P'],
                          target_name,
                          target_sequence)