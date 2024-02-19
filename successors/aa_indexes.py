__author__ = "Pavel Kohout xkohou15@stud.fit.vutbr.cz"
__date__ = "08-15-2022"
__description__ = "prepare AA indexes ranges"


import pickle
import numpy as np
import pandas as pd
import seaborn as sn
from scipy.spatial import distance_matrix

aa_idx1 = {
    "record":  "DAWD720101",
    "name": "size.pkl",
    "A": 2.5, "L": 5.5, "R": 7.5, "K": 7.0, "N": 5.0, "M": 6.0, "D": 2.5, "F": 6.5, "C": 3.0, "P": 5.5,
    "Q": 6.0, "S": 3.0, "E": 5.0, "T": 5.0, "G": 0.5, "W": 7.0, "H": 6.0, "Y": 7.0, "I": 5.5, "V": 5.0
}

aa_idx2 = {
    "record":  "FASG760101",
    "name": "molecular_weight.pkl",
    "A": 89.09, "L": 131.17, "R": 174.20, "K": 146.19, "N": 132.12, "M": 149.21, "D": 133.10, "F": 165.19, "C": 121.15, "P": 115.13,
    "Q": 146.15, "S": 105.09, "E": 147.13, "T": 119.12, "G": 75.07, "W": 204.24, "H": 155.16, "Y": 181.19, "I": 131.17, "V": 117.15
}

aa_idx3 = {
    "record":  "FASG760102",
    "name": "melting_point.pkl",
    "A": 297, "L": 337, "R": 238, "K": 224, "N": 236, "M": 283, "D": 270, "F": 284, "C": 178, "P": 222,
    "Q": 185, "S": 228, "E": 249, "T": 253, "G": 290, "W": 282, "H": 277, "Y": 344, "I": 284, "V": 293
}

aa_idx4 = {
    "record":  "FASG760102",
    "name": "waals_volume.pkl",
    "A": 1.0, "L": 4.0, "R": 6.13, "K": 4.77, "N": 2.95, "M": 4.43, "D": 2.78, "F": 5.89, "C": 2.43, "P": 2.72,
    "Q": 3.95, "S": 1.60, "E": 3.78, "T": 2.6, "G": 0.0, "W": 8.08, "H": 4.66, "Y": 6.47, "I": 4.0, "V": 3.0
}

aa_idx5 = {
    "record":  "GOLD730102",
    "name": "residue_volume.pkl",
    "A": 88.3, "L": 168.5, "R": 181.2, "K": 175.6, "N": 125.1, "M": 162.2, "D": 110.8, "F": 189.0, "C": 112.4, "P": 122.2,
    "Q": 148.7, "S": 88.7, "E": 140.5, "T": 118.2, "G": 60.0, "W": 227.0, "H": 152.6, "Y": 193.0, "I": 168.5, "V": 141.4
}

aa_idx6 = {
    "record":  "WOLR790101",
    "name": "hydrophobicity.pkl",
    "A": 1.12, "L": -2.55, "R": -2.55, "K": -0.8, "N": -0.83, "M": 0.55, "D": -0.83, "F": 0.67, "C": 0.59, "P": 0.54,
    "Q": -0.78, "S": -0.05, "E": -0.92, "T": -0.02, "G": 1.20, "W": -0.19, "H": -0.93, "Y": -0.23, "I": 1.16, "V": 1.13
}

aa_idx7 = {
    "record":  "KLEP840101",
    "name": "net_charge.pkl",
    "A": 0.0, "L": 0.0, "R": 1.0, "K": 1.0, "N": 0.0, "M": 0.0, "D": -1.0, "F": 0.0, "C": 0.0, "P": 0.0,
    "Q": 0.0, "S": 0.0, "E": -1.0, "T": 0.0, "G": 0.0, "W": 0.0, "H": 0.0, "Y": 0.0, "I": 0.0, "V": 0.0
}

aa_idx8 = {
    "record":  "BHAR880101",
    "name": "flexibility.pkl",
    "A": 0.357, "L": 0.365, "R": 0.529, "K": 0.466, "N": 0.463, "M": 0.295, "D": 0.511, "F": 0.314, "C": 0.346, "P": 0.509,
    "Q": 0.493, "S": 0.507, "E": 0.497, "T": 0.444, "G": 0.544, "W": 0.305, "H": 0.323, "Y": 0.42, "I": 0.462, "V": 0.386
}

aa_idx9 = {
    "record":  "BULH740101",
    "name": "free_energy_surface.pkl",
    "A": -0.2, "L": -2.46, "R": -0.12, "K": -0.35, "N": 0.08, "M": -1.47, "D": -0.2, "F": -2.33, "C": -0.45, "P": -0.98,
    "Q": 0.16, "S": -0.39, "E": -0.3, "T": -0.52, "G": 0.0, "W": -2.01, "H": -0.12, "Y": -2.24, "I": -2.26, "V": -1.56
}

aa_idx10 = {
    "record":  "FAUJ880108",
    "name": "localized_electrical_effect.pkl",
    "A": -0.01, "L": -0.01, "R": 0.04, "K": 0.00, "N": 0.06, "M": 0.04, "D": 0.15, "F": 0.03, "C": 0.12, "P": 0.0,
    "Q": 0.05, "S": 0.11, "E": 0.07, "T": 0.04, "G": 0.0, "W": 0.0, "H": 0.08, "Y": 0.03, "I": -0.01, "V": 0.01
}

aa_idx11 = {
    "record":  "ZIMJ680103",
    "name": "polarity.pkl",
    "A": 0.00, "L": 0.13, "R": 52.0, "K": 49.5, "N": 3.38, "M": 1.43, "D": 49.7, "F": 0.35, "C": 1.48, "P": 1.58,
    "Q": 3.53, "S": 1.67, "E": 49.9, "T": 1.66, "G": 0.00, "W": 2.1, "H": 51.6, "Y": 1.61, "I": 0.13, "V": 0.13
}

aa_idx12 = {
    "record":  "ZIMJ680104",
    "name": "isoelectric.pkl",
    "A": 6.0, "L": 5.98, "R": 10.76, "K": 9.74, "N": 5.41, "M": 5.74, "D": 2.77, "F": 5.48, "C": 5.05, "P": 6.3,
    "Q": 5.65, "S": 5.68, "E": 3.22, "T": 5.66, "G": 5.97, "W": 5.89, "H": 7.59, "Y": 5.66, "I": 6.02, "V": 5.96
}

aa_idxs = [aa_idx2, aa_idx3, aa_idx5, aa_idx6, aa_idx8, aa_idx9, aa_idx10, aa_idx11, aa_idx12]

# for each index create matrix of distances and percentage similarities
alphabet = ["A", "L", "R", "K", "N", "M", "D", "F", "C", "P", "Q", "S", "E", "T", "G", "W", "H", "Y", "I", "V"]
alphabet_idx = {}
for i, a in enumerate(alphabet):
    alphabet_idx[a] = i
with open("./aa_index_confidence/aa_idx.pkl", "wb") as f:
    pickle.dump(alphabet_idx, f)


for aa_idx in aa_idxs:
    aa_vals = [[aa_idx[k]] for k in alphabet]
    matrix = distance_matrix(aa_vals, aa_vals)
    matrix = 1.0 - (matrix / np.max(matrix))

    with open("./aa_index_confidence/" + aa_idx["name"], "wb") as f:
        pickle.dump(matrix, f)
    print("=" * 80)
    print(aa_idx["name"])
    print()

file_names = [d["name"] for d in aa_idxs]
with open("./aa_index_confidence/file_names.pkl", "wb") as f:
    pickle.dump(file_names, f)

indices = [d["record"] for d in aa_idxs]
with open("./aa_index_confidence/aa_indice_names.pkl", "wb") as f:
    pickle.dump(indices, f)

data = dict()
for ind in aa_idxs:
    # name = ind['name'].split(".")[0]
    name = ind['record']
    data[name] = [ind[v] for v in alphabet]

df = pd.DataFrame(data)
corr_matrix = df.corr()
svm = sn.heatmap(corr_matrix, annot=True)
svm.figure.tight_layout()
figure = svm.get_figure()
figure.savefig('./aa_index_confidence/corr_matrix_ind.png', dpi=400)
