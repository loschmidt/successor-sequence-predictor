__author__ = "Pavel Kohout <xkohou15@vutbr.cz>, Milos Musil <imusil@vut.cz>"
__date__ = "2024/02/16"
__description__ = "Predict successor amino acids for all features (AA indices you have selected)"

import glob
import itertools
import os
import pickle
from typing import Tuple
from successors.internal_types import *

import numpy as np
from Bio import Phylo
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

from successors.parser_handler import RunSetup
from successors.utils import load_msa

target_vec = ""


def log_prediction(run_config: RunSetup, trees: int) -> str:
    msg = "#" * 80 + "\n"
    msg += f'   Predicting successors of {run_config.protein_name} with {trees} trees\n\n'
    msg += '       predicting ' + ('extant sequence of ' if run_config.validation else 'successor sequence of ')
    msg += f' {run_config.protein_name} with tree identifier {run_config.query}\n'
    msg += '       ' + 'tree lineages with more than three transitions between AAs will be used!' \
        if run_config.transition else 'all tree lineages are used!'
    msg += "\n"
    return msg


def linage_trend(feature_vals):
    """ Create sequencetiality measurement for lineage, return sequenciality, sequenciality2, fluctuation"""
    # group by following same AAs, e.g. AABBBCDD -> ABCD
    grouped = [feature_vals[0]] + [val for i, val in enumerate(feature_vals[1:]) if val != feature_vals[i]]
    fluctuation = 100 * len(set(grouped)) / len(grouped)  # originality of AA in lineage

    # for every pair in grouped detect raising/descending/neutral trend +1/-1/+0
    sequential_trend_sum, pairs = 0, list(itertools.combinations(grouped, 2))
    for a, b in pairs:
        if a > b:
            sequential_trend_sum -= 1
            continue
        if a < b:
            sequential_trend_sum += 1
            continue

    if len(pairs) == 0:
        return 0, 0, 0
    sequential_trend = (100 * sequential_trend_sum) / len(pairs)

    # check whether there are break in trend for last transitions
    if len(grouped) < 3:
        sequential_trend_break = -1  # nothing to detect
    else:
        a, b, c = grouped[-3:]
        sequential_trend_break = ((a < b) and (b > c)) or ((a > b) and (b < c))
    return sequential_trend, sequential_trend_break, fluctuation


def predict_tree_successor(lineage_sequences: LINEAGE_SEQUENCES, features: AA_INDICES, tree_dir: str, run: RunSetup):
    """
    Predict tree successor with selected features and collect important statistics.
    :param lineage_sequences: Parsed sequences from the tree in a path from query to the root, list per position
    :param features: Features dictionary
    :param tree_dir: Tree directory to generate tmp file to
    :param run: RunSetup object with experiment configuration
    :return:
    """
    global target_vec

    msa_width = len(lineage_sequences)
    ft_list = list()
    for i in range(msa_width):
        ft_list.append(dict())

    scores = {"r2": {f: [None for _ in range(msa_width)] for f in features},
              "slope": {f: [None for _ in range(msa_width)] for f in features},
              "intercept": {f: [None for _ in range(msa_width)] for f in features},
              "values": {f: [None for _ in range(msa_width)] for f in features},
              "aminos": {f: [None for _ in range(msa_width)] for f in features},
              "pred_vals": {f: [None for _ in range(msa_width)] for f in features},
              "target_vec": {f: [None for _ in range(msa_width)] for f in features},
              "sequentiality": {f: [None for _ in range(msa_width)] for f in features},
              "fluctuation": {f: [None for _ in range(msa_width)] for f in features},
              "sequentiality2": {f: [None for _ in range(msa_width)] for f in features}}

    positions_to_keep, positions_detected = [], False
    aa_on_the_last = []
    for feature in features:
        for i in range(msa_width):
            act = list()
            predict_from_feat = list()

            # get columns for prediction and translate AA to feature values
            transitions_list = list()
            for vc in lineage_sequences[i]:
                if vc == '-':  # sequence has the gap on this MSA column skip it
                    predict_from_feat.append("-")
                    continue
                else:
                    act.append(float(features[feature][vc]))
                    transitions_list.append(vc)
                    predict_from_feat.append(float(features[feature][vc]))

            # keep tract of AA in the last ancestor
            aa_on_the_last.append(lineage_sequences[i][-1])
            # check the youngest ancestor, keep a gap
            if predict_from_feat[-1] == "-":
                ft_list[i][feature] = "-"
                continue

            if not positions_detected:
                positions_to_keep.append(i)
            # aggregate transitions and allow predictions for 3+ otherwise keep youngest
            transitions_cnt, last_transition = 0, transitions_list[0]
            for tr in transitions_list[1:]:
                if tr != last_transition:
                    last_transition = tr
                    transitions_cnt += 1

            # do prediction unless less than 3 transitions occurred for position and transition check is on
            if transitions_cnt < 3 and run.transition:
                ft_list[i][feature] = predict_from_feat[-1]
                continue

            # do prediction only if features in position is more than half of the tree depth to the query
            if len(act) * 2 >= len(lineage_sequences[i]):
                x = list()
                y_w = list()
                vals_inds = list()

                # weights for prediction, set weight linearly to the distance from the root of the tree
                for j in range(len(act)):
                    for weight in range(j + 1):
                        x.append(len(y_w))
                        y_w.append(act[j])
                        vals_inds.append(x[-1] - 1)
                x = np.array(x).reshape((-1, 1))
                y = np.array(y_w)
                model = LinearRegression().fit(x, y)
                y_pred_all = model.predict(x)

                # predict feature value for last, previous x-axis values are from 0 so now plus one (len does it)
                y_pred = model.predict(np.array(len(y_w)).reshape(-1, 1))
                ft_list[i][feature] = y_pred

                scores['r2'][feature][i] = r2_score(y[vals_inds], y_pred_all[vals_inds])
                scores['slope'][feature][i] = model.coef_[0]
                scores['intercept'][feature][i] = model.intercept_
                scores['values'][feature][i] = y
                scores['aminos'][feature][i] = transitions_list
                scores['pred_vals'][feature][i] = [y_pred]
                scores['target_vec'][feature][i] = [float(features[feature][target_vec[i]]) if target_vec[i] != '-'
                                                    else 0.0, target_vec[i]]
                trend, trend_break, fluctuation = linage_trend(act)
                scores['sequentiality'][feature][i] = trend
                scores['sequentiality2'][feature][i] = trend_break
                scores['fluctuation'][feature][i] = fluctuation
            else:
                # otherwise, keep what is in youngest (the value of current feature)
                ft_list[i][feature] = predict_from_feat[-1]
        positions_detected = True  # detect non gap position just for the first feature, for all others it is the same

    # make prediction base on individual AA indices
    indices_pred = {k: ["", "", ""] for k in features if
                    k != "gap"}  # 1st 2nd 3rd probable

    aa_list = ["A", "L", "R", "K", "N", "M", "D", "F", "C", "P", "Q", "S", "E", "T", "G", "W", "H", "Y", "I", "V"]

    for i in range(msa_width):
        # feature wise aa selection
        for feature in ft_list[i]:
            keep_vals_aa = list()

            if str(ft_list[i][feature]) == "-":
                continue

            for aa in aa_list:  # get AA from given feature with the closest difference from prediction
                aa_dev = (float(ft_list[i][feature]) - float(features[feature][aa])) * (
                        float(ft_list[i][feature]) - float(features[feature][aa]))

                keep_vals_aa.append((aa, aa_dev))
            keep_vals_aa.sort(key=lambda xx: xx[1])
            if scores['pred_vals'][feature][i] is not None:
                scores['pred_vals'][feature][i].append(keep_vals_aa[0][0])
            indices_pred[feature][0] += keep_vals_aa[0][0]
            indices_pred[feature][1] += keep_vals_aa[1][0]
            indices_pred[feature][2] += keep_vals_aa[2][0]

            # For some indices, there is more AA sharing same indices values.
            # Therefore, check if the original AA symbol is along those with the best prediction.
            # If so, keep the original AA as the prediction to remove AA alphabet order influence and reduce
            # false substitutions in designs
            the_best_aa, lowest_dev = [], keep_vals_aa[0][1]
            for b in keep_vals_aa:
                if b[1] <= lowest_dev:
                    the_best_aa.append(b[0])
            if aa_on_the_last[i] in the_best_aa:
                indices_pred[feature][0] = indices_pred[feature][0][:-1] + aa_on_the_last[i]  # replace that last value
                if scores['pred_vals'][feature][i] is not None:
                    scores['pred_vals'][feature][i][-1] = aa_on_the_last[i]  # replace that last value

    # create a file for individual candidates 3 best candidates
    for pred_i in range(3):
        prediction_file_path = os.path.join(tree_dir, "prediction_indices_{}.fasta".format(pred_i))
        ot = open(prediction_file_path, "w")
        for k, seq in indices_pred.items():
            ot.write(">{}\n".format(k))
            ot.write(seq[pred_i])
            ot.write("\n")
        ot.close()

    # remove all positions with gap as results
    new_scores = {}
    for m in scores:
        new_scores[m] = {}
        for feature, vals in scores[m].items():
            new_scores[m][feature] = [vals[i_good] for i_good in positions_to_keep]
    scores = new_scores
    # store slopes and r2
    for m in scores:
        for feature, vals in scores[m].items():
            feature_file_path = os.path.join(tree_dir, "{}_{}.pkl".format(m, feature))
            with open(feature_file_path, "wb") as f:
                pickle.dump(vals, f)


def tree_wise_predictions(features: AA_INDICES, trees: TREE_PATHS, run: RunSetup) -> str:
    """
    Predict tree ancestor for each tree in an input dataset
    :param features: normalized feature map
    :param trees: paths with given trees in input dataset
    :param run: configuration of the run
    :return:
    """
    msg = "Prediction done for:"

    for i, tree_dir_path in enumerate(trees):
        print("Predicting " + tree_dir_path)
        # create temporary experiment run location
        tree_tmp_dir = os.path.join(run.tmp_files, f"tree{i}")
        os.makedirs(tree_tmp_dir, exist_ok=True)

        try:
            lineage_sequences = get_tree_lineage_sequences(tree_dir_path, run)
        except FileNotFoundError as e:
            print(f"Not found file in tree input dir: {tree_dir_path}!!!")
            print(e)
            exit(1)
        except Exception as e:
            print("Clade or different error while parsing sequences in the phylogeny tree lineage")
            print(e)
            exit(1)

        predict_tree_successor(lineage_sequences, features, tree_tmp_dir, run)

        # Update log message
        msg += f" {i}-{tree_dir_path},"
        if (i + 1) % 6 == 0:
            msg += "\n             "

    return msg


def get_tree_lineage_sequences(dataset_tree_path: str, run: RunSetup) -> LINEAGE_SEQUENCES:
    # sequences = dict()
    # length = 0
    # dt = open(dataset_tree_path + "/msa.fasta")
    # name = ""
    # act = ""
    # for line in dt.readlines():
    #     line = line.strip()
    #     if line[0] == '>':
    #         if name != "":
    #             sequences[name] = act
    #             if length == 0:
    #                 length = len(act)
    #         name = line
    #         act = ""
    #     else:
    #         act += line
    # sequences[name] = act

    sequences, length = load_msa(os.path.join(dataset_tree_path, "msa.fasta"))

    tree = Phylo.read(dataset_tree_path + "/ancestralTree.tree", "newick")

    global target_vec

    query_name = run.query
    nodes = list()
    for clade in tree.get_terminals():
        path = tree.get_path(clade)
        if query_name.lower() in str(path[-1]).lower():
            for cl in path:
                nodes.append(">ancestral_" + str(cl.confidence))
            nodes[-1] = ">" + path[-1].name
            target_vec = sequences[nodes[-1]]
            if run.validation:
                nodes = nodes[:-1]  # hide leave node as target TODO check it twice
            break
    # get a list of AA symbols from root to the query leaf for each column in MSA
    vectors = list()
    for i in range(length):
        vectors.append(list())
        for node in nodes:
            vectors[i].append(sequences[node][i])

    return vectors


def predict(run_config: RunSetup) -> None:
    """
    Predict successor sequences from run configuration passed via parameters
    :param run_config:
    :return:
    """
    input_pattern = os.path.join(run_config.input, '*')
    trees = glob.glob(input_pattern)

    log_msg = log_prediction(run_config, len(trees))
    print(log_msg)
    log_file_path = os.path.join(run_config.logs, "log_file.txt")
    log_file = open(log_file_path, "w")
    log_file.write(log_msg)

    # load AA indices features
    features, log_msg = get_indices_features(run_config)
    print(log_msg)
    log_file.write(log_msg)

    # make predictions of successor with given trees
    log_msg = tree_wise_predictions(features, trees, run_config)
    print(log_msg)
    log_file.write(log_msg)


def get_indices_features(run_config: RunSetup) -> Tuple[AA_INDICES, str]:
    """
    Get normalized values of features used for the regression over trees
    :param run_config:
    :return: Dictionary with parsed normalized indices values, log message
    """
    aaindex_file = open(os.path.join('successors', 'indices', 'aaindex.csv'))

    used_indices = []
    unused_indices = []

    features = dict()
    first = True
    for line in aaindex_file:
        line = line.strip()
        if first:
            aa = line.split(",")[1:]
            first = False
        else:
            parts = line.split(",")
            data = parts[1:]
            name = parts[0]
            if name not in run_config.indices:
                unused_indices.append(name)
                continue
            used_indices.append(name)
            act = dict()
            for i in range(len(data)):
                act[aa[i]] = float(data[i])
            features[name] = act

    # normalize values of features
    for feature in features:
        min_f = None
        max_f = None
        for ft in features[feature]:
            if min_f == None or float(features[feature][ft]) < min_f:
                min_f = float(features[feature][ft])
            if max_f == None or float(features[feature][ft]) > max_f:
                max_f = float(features[feature][ft])
        for ft in features[feature]:
            act = float(features[feature][ft])
            norm = (act - min_f) / (max_f - min_f)
            features[feature][ft] = norm

    log_msg = "#" * 80 + "\n"
    log_msg += "   Getting features for....\n"
    log_msg += "     " + ",".join(used_indices) + "\n"
    log_msg += "   These features are available but were not used\n"
    log_msg += "     " + ",".join(unused_indices) + "\n\n"

    return features, log_msg
