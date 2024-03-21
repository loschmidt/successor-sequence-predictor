__author__ = "Pavel Kohout xkohou15@vutbr.cz"
__date__ = "02-19-2024"
__description__ = "3rd level of evaluation"

import os
import pickle

import pandas as pd

from successors.parser_handler import RunSetup

CONS_SCORE = None
CONS_POS = 0
CONS_OFF = 0


def load_conservation_from_file(file_path):
    first_line = True
    positions, seqs, scores = [], [], []
    with open(file_path) as f:
        for line in f:
            if first_line:
                first_line = False
                continue
            line = line.split()
            # we want position, sequence, score
            pos, seq, score = line[0], line[1], line[4]
            positions.append(pos)
            seqs.append(seq)
            scores.append(score)
    return positions, seqs, scores


def add_table_head(final_dict, title, columns):
    """ Add header into generated table with title 3A for instance """
    final_dict[0].extend(["", title, columns[0]])
    for k in list(final_dict.keys())[1:]:
        final_dict[k].extend(["", "", columns[k]])
    return final_dict


def add_table_records(final_dict, mutations, level):
    """ Mutations per positions sorted by count of them """
    global CONS_SCORE
    global CONS_OFF
    global CONS_POS

    for muts in mutations:
        muts = sorted(muts, key=lambda x: (x[2], x[3]), reverse=True)
        for j, mut in enumerate(muts):
            sub, tree, seq, fluc, br, N, G, AAC, msa_i, indice = mut
            if level == "a":
                if j == 0:
                    final_dict[0].extend([sub[:-1]])  # substitution
                    final_dict[1].extend([sub[-1]])  # change to
                else:
                    final_dict[0].extend([""])  # substitution
                    final_dict[1].extend([""])  # change to
            if level == "b":
                if j == 0:
                    final_dict[0].extend([sub[:-1]])  # substitution
                else:
                    final_dict[0].extend([""])
                final_dict[1].extend([sub[-1]])  # change to
            if level == "c":
                final_dict[0].extend([sub[:-1]])  # substitution
                final_dict[1].extend([sub[-1]])  # change to
            # actually append record to column
            if int(CONS_POS[-1]) > msa_i >= CONS_OFF:
                cons_i = msa_i - CONS_OFF
                cons = CONS_SCORE[cons_i]
            else:
                cons = ""
            for mut_i, item in enumerate([indice, tree, cons, seq, fluc, br, N, G, AAC]):
                final_dict[2 + mut_i].extend([item])
    return final_dict


def a3_table(ind_stats, final_dict, columns):
    """Apply a3 level search over indices same position and same AA"""
    final_dict = add_table_head(final_dict, "3A", columns)

    over_ind = {}
    # select over indices without break
    for name, stats in ind_stats.items():
        for stat in stats:
            if len(stat) == 0: continue
            if stat[4] == 1 or stat[4] == -1:
                continue
            sub, tree, seq, fluc, br, N, G, AAC, cons = stat
            if not (sub in over_ind.keys()):
                over_ind[sub] = []
            add_record = stat + [name]  # add indice name to record
            over_ind[sub].append(add_record)

    # pick the best tree from each indice
    one_tree_per_indice = {}
    for substitutions in over_ind.values():
        indices = list(set([rec[-1] for rec in substitutions]))
        substitution = substitutions[0][0]  # get identifier of mutations
        one_tree_per_indice[substitution] = []
        for indice in indices:
            indice_subs = [rec for rec in substitutions if rec[-1] == indice]
            best_per_ind = sorted(indice_subs, key=lambda x: (x[2], x[3]))[-1]
            one_tree_per_indice[substitution].append(best_per_ind)
    over_ind = one_tree_per_indice

    # sort by predictions counts
    sort_over_ind = sorted(over_ind, key=lambda k: len(over_ind[k]), reverse=True)
    sort_over_ind = [over_ind[k] for k in sort_over_ind]
    final_dict = add_table_records(final_dict, sort_over_ind, "a")
    return final_dict


def b3_table(ind_stats, final_dict, columns):
    """Apply b3 level search over indices same position and different AA"""
    final_dict = add_table_head(final_dict, "3B", columns)

    over_ind = {}
    # select over indices without break
    for name, stats in ind_stats.items():
        for stat in stats:
            if len(stat) == 0: continue
            if stat[4] == 1 or stat[4] == -1:
                continue
            sub, tree, seq, fluc, br, N, G, AAC, cons = stat
            sub = sub[:-1]
            if not (sub in over_ind.keys()):
                over_ind[sub] = []
            add_record = stat + [name]  # add indice name to record
            over_ind[sub].append(add_record)

    # pick the best tree from each indice
    one_tree_per_indice = {}
    for substitutions in over_ind.values():
        indices = list(set([rec[-1] for rec in substitutions]))
        substitution = substitutions[0][0]  # get identifier of mutations
        one_tree_per_indice[substitution] = []
        for indice in indices:
            indice_subs = [rec for rec in substitutions if rec[-1] == indice]
            best_per_ind = sorted(indice_subs, key=lambda x: (x[2], x[3]))[-1]
            one_tree_per_indice[substitution].append(best_per_ind)
    over_ind = one_tree_per_indice

    # sort by predictions counts
    sort_over_ind = sorted(over_ind, key=lambda k: len(over_ind[k]), reverse=True)
    sort_over_ind = [over_ind[k] for k in sort_over_ind]
    final_dict = add_table_records(final_dict, sort_over_ind, "b")
    return final_dict


def predict_level3(run: RunSetup):
    """
    Predict level three sequence
    :param run:
    :return:
    """
    columns = ["Residue number", "Mutated to", "According to Index", "Tree", "Consurf", "Sequentiality", "Fluctuation",
               "Break trend", "N", "G", "AA CNT"]

    level3_dir = os.path.join(run.results, "level3")  # f"outputs/{protein}/level3/"
    os.makedirs(level3_dir, exist_ok=True)

    level2_dir = os.path.join(run.results, "level2")
    CONFIDENCE_LEVEL = run.confidence_level

    protein = run.protein_name
    # for protein in proteins:

    print(level3_dir, "CONFIDENCE LEVEL {}".format(CONFIDENCE_LEVEL))

    gen_dir_level2 = f"outputs/{protein}/level2/script_tmp/"
    indices_path = os.path.join(level2_dir, f"indices_transitions{CONFIDENCE_LEVEL}.pkl")
    # with open(gen_dir_level2 + "indices_transitions{}.pkl".format(CONFIDENCE_LEVEL), "rb") as f:
    with open(indices_path, "rb") as f:
        ind_stats = pickle.load(f)

    global CONS_POS
    global CONS_SCORE

    CONS_POS = [0 for _ in range(1000)]
    CONS_SCORE = CONS_POS

    # init final dict to csv columns
    final_dict = {i: [] for i in range(len(columns))}
    final_dict = a3_table(ind_stats, final_dict, columns)

    df = pd.DataFrame.from_dict(final_dict)
    df.to_csv(os.path.join(level3_dir, "a3_{}_metrics_{}.csv".format(protein, CONFIDENCE_LEVEL)))

    # init final dict to csv columns
    final_dict = {i: [] for i in range(len(columns))}
    final_dict = b3_table(ind_stats, final_dict, columns)

    df = pd.DataFrame.from_dict(final_dict)
    df.to_csv(os.path.join(level3_dir, "b3_{}_metrics_{}.csv".format(protein, CONFIDENCE_LEVEL)))
