__author__ = "Pavel Kohout xkohou15@vutbr.cz"
__date__ = "20/02/2024"
__description__ = "Keep often used logic together"

import os

from successors.internal_types import MSA_GAP_LEN, MSA


def load_msa(msa_file_path: str) -> MSA_GAP_LEN:
    sequences = dict()
    msa_file = open(msa_file_path, "r")
    length = 0
    name = ""
    act = ""
    for line in msa_file.readlines():
        line = line.strip()
        if line[0] == '>':
            if name != "":
                sequences[name] = act
                if length == 0:
                    length = len(act)
            name = line
            act = ""
        else:
            act += line
    sequences[name] = act
    return sequences, length


def store_msa(msa_file_path: str, sequences: MSA):
    with open(msa_file_path, "w") as f:
        for n, s in sequences.items():
            f.write(f">{n}\n{s}\n")


def clustalo_sequences(to_align_path: str, msa_path: str, sequences: MSA) -> MSA_GAP_LEN:
    store_msa(to_align_path, sequences)
    os.system(f"clustalo -i {to_align_path} -o {msa_path} --force")
    return load_msa(msa_path)
