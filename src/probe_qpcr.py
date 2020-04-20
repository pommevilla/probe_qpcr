def read_primers(file_name: str) -> dict:
    """
    Reads primers from file.
    Returns a dictionary of dictionaries of lists.
    """
    primer_dict = {}
    with open(file_name) as fin:
        name, seq, direction = None, "", None
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    primer_dict[name][direction] = seq
                *name, direction = line.split(".")
                # TODO: Is there a better way to handle primer names?
                name[0] = name[0][1:]
                name[-1] = name[-1][:-2]
                name = '.'.join(name)
                seq = ""
                if not primer_dict.get(name):
                    primer_dict[name] = {'F': '', 'R': '', 'P': ''}
            else:
                seq = replace_ambiguity_bases(line)
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


def read_sequences(file_name: str) -> iter:
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
                name, seq = line[1:], []
            else:
                seq.append(line)
        if name:
            yield name, ''.join(seq)


# TODO: Change method to accept a whole dictionary entry
def get_qpcr_hits(primer_name: str, forward_primer: str, reverse_primer: str, probe: str,
                  target_name: str, target_sequence: str) -> str:
    """
    Returns the TaqMan product.
    """
    import re
    primer_pattern = "({}).*({}).*({})".format(forward_primer, probe, reverse_complement(reverse_primer))
    for match in [match for match in re.finditer(primer_pattern, target_sequence)]:
        product = target_sequence[match.start():match.end()]
        print("{}\t{}\t{}\t{}\t{}\t{}".format(primer_name, target_name, match.start(), match.end(),
                                              match.end() - match.start(), product))


def replace_ambiguity_bases(sequence: str) -> str:
    codes = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "S": "[GC]", "B": "[CGT]", "Y": "[CT]", "W": "[AT]",
             "D": "[AGT]", "K": "[GT]", "N": "[ACGT]", "H": "[ACT]", "M": "[AC]", "V": "[ACG]", "X": "[ACGT]"}
    regexlist = []
    for base in sequence:
        try:
            regexlist.append(codes[base])
        except KeyError:
            print(f"Base {base} not found.")
            raise

    return ''.join(regexlist)


def main():
    primer_file, target_file = sys.argv[1:]
    primer_dict = read_primers(primer_file)
    for target_name, target_sequence in read_sequences(target_file):
        for primer_set in primer_dict:
            # TODO: change method signature
            get_qpcr_hits(primer_set,
                          primer_dict[primer_set]['F'],
                          primer_dict[primer_set]['R'],
                          primer_dict[primer_set]['P'],
                          target_name,
                          target_sequence)


if __name__ == "__main__":
    # TODO: Change to sys.exit() syntax
    import sys
    sys.exit(main())
