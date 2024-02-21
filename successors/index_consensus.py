__author__ = "Pavel Kohout xkohou15@vutbr.cz"
__date__ = "2024/02/19"
__description__ = "For every AA indice make an average and compute WT similarity over all trees"

import glob
import math
import os
import pickle
import statistics

import pandas as pd

from successors.parser_handler import RunSetup
from successors.utils import load_msa, clustalo_sequences


def create_stats_with_metric(vals_over_trees, trees_msa, msa_mapping, msa_width):
    """
    Calculates statistic distribution of values in metrics
    :param vals_over_trees: array of arrays with values per sequence position (None included)
    :param trees_msa: indice alignment over trees
    :param msa_mapping: dict of keys seq_name and mapping of sequence residues to msa positions
    :param msa_width: width of msa
    :return: 3 collections,
            metric_stats - list of list having 5 elements with r2 mean and so on
            vals_per_tree_per_col - list of list of list with 2 elem [slope/r2 val, aa_char]
            voting_aa_over_prediction - list for columns, None = gap col, [] no predictions for col,
                                        [elems_aa] aa chars predicted by regression
    """
    metric_stats = list()

    # reduce None positions in
    msa_seq_names = list(msa_mapping.keys())  # first one is a query
    # print(len(vals_over_trees), len(msa_seq_names))
    mapping_folder = list()
    for i, name_al_seq in enumerate(msa_seq_names):
        mapping_folder.append(msa_mapping[name_al_seq])
    vals_per_tree_per_col = list()
    voting_aa_over_prediction = list()  # consensus based on predictions only
    for pos_i in range(msa_width):
        tree_per_col = list()
        prediction_aa_list = list()
        for f_i, v in enumerate(mapping_folder):
            aa_in_tree = trees_msa[msa_seq_names[f_i]][pos_i]
            if pos_i in v:
                val = vals_over_trees[f_i][v.index(pos_i)]
                if val:
                    tree_per_col.append([val, aa_in_tree])
                    prediction_aa_list.append(aa_in_tree)
                else:
                    tree_per_col.append([math.nan, aa_in_tree])
            else:
                tree_per_col.append([math.nan, aa_in_tree])
                prediction_aa_list.append(None)
        vals_per_tree_per_col.append(tree_per_col)
        voting_aa_over_prediction.append(prediction_aa_list)

        vals_with_none = [vals_over_trees[f_i][v.index(pos_i)] for f_i, v in enumerate(mapping_folder) if pos_i in v]
        vals = [v for v in vals_with_none if v]
        if len(vals) == 0:
            metric_stats.append([0.0, 0.0, 0.0, 0.0, 0])
            continue

        cnt = len(vals)
        met_max, met_min = max(vals), min(vals)
        mean = statistics.mean(vals)
        if cnt >= 2:
            met_stddev = statistics.stdev(vals)
        else:
            met_stddev = 0.0
        metric_stats.append([mean, met_stddev, met_max, met_min, cnt])

    return metric_stats, vals_per_tree_per_col, voting_aa_over_prediction


def voting_by_prediction(predictions):
    """voted prediction for AA"""
    predictions = [p for p in predictions if p is not None]  # remove Nones
    if predictions is None or predictions == []:
        return '', 0
    freq = {aa: 0 for aa in set(predictions)}
    for prediction in predictions:
        freq[prediction] += 1
    ordered_freq = {k: v for k, v in sorted(freq.items(), key=lambda item: item[1])}
    the_best = list(ordered_freq.items())[-1]
    return the_best


def sequence_trend_metrics(metrics, feature, tree_map_to_msa, indice_msa_len):
    seq_data = metrics["sequence_trend"][feature]
    break_trend = metrics["sequence_trend_break"][feature]
    fluctuation = metrics["fluctuations"][feature]
    ys = metrics["values"][feature]

    templates = []
    for msa_pos in range(indice_msa_len):
        # S,F,B,P,G, plot
        template = [[0, 0, -1, 0, 0, -1] for _ in range(TREE_CNT)]
        for tree_i, tree_msa_pos in enumerate(tree_map_to_msa):
            if msa_pos not in tree_msa_pos:
                continue
            if seq_data[tree_i][tree_msa_pos.index(msa_pos)] is None:  # was there the prediction?
                continue  # if not continue
            # pos_over_sequential.append(seq_data[tree_i][tree_msa_pos.index(msa_pos)])
            # pos_over_break.append(break_trend[tree_i][tree_msa_pos.index(msa_pos)])
            # pos_over_fluct.append(fluctuation[tree_i][tree_msa_pos.index(msa_pos)])
            template[tree_i][0] = seq_data[tree_i][tree_msa_pos.index(msa_pos)]
            template[tree_i][1] = fluctuation[tree_i][tree_msa_pos.index(msa_pos)]
            template[tree_i][2] = break_trend[tree_i][tree_msa_pos.index(msa_pos)]

            feature_vals = ys[tree_i][tree_msa_pos.index(msa_pos)]
            grouped = [feature_vals[0]] + [val for i, val in enumerate(feature_vals[1:]) if val != feature_vals[i]]

            template[tree_i][3] = len(ys[tree_i][tree_msa_pos.index(msa_pos)])
            template[tree_i][4] = len(grouped)
            template[tree_i][5] = tree_msa_pos.index(msa_pos)
        templates.append(template)
    return templates


def create_stats_with_metric2(trees_metric, trees_msa, tree_map_to_msa, msa_width):
    """
        Calculates statistic distribution of values in metrics
        :param trees_metric: array of arrays with values per sequence position (None included)
        :param trees_msa: indice alignment over trees
        :param tree_map_to_msa: list of msa position to sequence position translation
        :param msa_width: width of msa
        :return: 3 collections,
                metric_stats - list of list having 5 elements with r2 mean and so on
                vals_per_tree_per_col - list of list of list with 2 elem [slope/r2 val, aa_char]
                voting_aa_over_prediction - list for columns, None = gap col, [] no predictions for col,
                                            [elems_aa] aa chars predicted by regression
    """
    metric_stats = list()
    vals_per_tree_per_col = list()
    voting_aa_over_prediction = list()  # consensus based on predictions only

    # convert MSA dict to array
    msa_index_list = list()
    for _, sequence in trees_msa.items():
        msa_index_list.append(sequence)

    # over position of MSA
    for msa_pos in range(msa_width):
        trees_per_col, prediction_aa_list, column_values = [], [], []
        there_was_prediction = False
        # prediction part
        for tree_i, tree_msa_pos in enumerate(tree_map_to_msa):
            aa_in_tree = msa_index_list[tree_i][msa_pos]
            if msa_pos not in tree_msa_pos:
                trees_per_col.append([math.nan, aa_in_tree])
                prediction_aa_list.append(None)
                continue
            tree_prediction_value = trees_metric[tree_i][tree_msa_pos.index(msa_pos)]
            if tree_prediction_value is None:  # was there the prediction?
                trees_per_col.append([math.nan, aa_in_tree])
                prediction_aa_list.append(None)
                continue  # if not continue

            there_was_prediction = True
            column_values.append(tree_prediction_value)
            trees_per_col.append([tree_prediction_value, aa_in_tree])
            prediction_aa_list.append(aa_in_tree)

        vals_per_tree_per_col.append(trees_per_col)
        voting_aa_over_prediction.append(prediction_aa_list)

        if not there_was_prediction:
            metric_stats.append([0.0, 0.0, 0.0, 0.0, 0])
            continue

        cnt = len(column_values)
        met_max, met_min = max(column_values), min(column_values)
        mean = statistics.mean(column_values)
        if cnt >= 2:
            met_stddev = statistics.stdev(column_values)
        else:
            met_stddev = 0.0
        metric_stats.append([mean, met_stddev, met_max, met_min, cnt])

    return metric_stats, vals_per_tree_per_col, voting_aa_over_prediction


def detect_best_metric(trends, trees_msa, aligned_wt, dict_line_cnt, protein, conservations):
    def calc_aa_in_lineage(n):
        aa_cnt, cumsum = 0, 0
        while True:
            cumsum += aa_cnt
            if cumsum == n:
                return aa_cnt
            aa_cnt += 1

    # conservation
    if not conservations:
        conservations = [0 for _ in range(1000)]  # fill it to make code clear
    else:  # apply offset
        off_cons = {"P0ACI0": 2, "P0A8U6": 1}[protein]
        conservations = ['' for _ in range(off_cons)] + conservations

    # convert MSA dict to array
    msa_index_list = list()
    for _, sequence in trees_msa.items():
        msa_index_list.append(sequence)

    best_dict = {"Individul trees": ["" for _ in range(dict_line_cnt)],
                 "Order": [], "Substitution": [], "Tree": [], "Sequentiality": [], "Fluctuation": [], "Break trend": [],
                 "N": [], "G": [], "AA CNT": [], "Conservation Score": []}
    p0a8u6 = {"P0A8U6": ["" for _ in range(dict_line_cnt)],
              "Order_2": [], "Substitution_2": [], "Tree_2": [], "Sequentiality_2": [], "Fluctuation_2": [],
              "Break trend_2": [],
              "N_2": [], "G_2": [], "AA CNT_2": [], "Conservation Score": []}
    # S,F,B,P,G
    to_sort = list()
    regardless_break = list()
    indice_stats = list()
    to_sort_pos = list()
    for msa_i, trees_data in enumerate(trends):
        for tree_i, data in enumerate(trees_data):
            # substitutions, tree
            if aligned_wt[msa_i] == msa_index_list[tree_i][msa_i]:
                continue  # that is not a substitution
            def_list = [f"{aligned_wt[msa_i]}{msa_i}{msa_index_list[tree_i][msa_i]}", tree_i]
            def_list.extend(data)
            def_list.extend([msa_i])
            indice_stats.append(def_list)
            # positions targets
            if protein == "P0A8U6":
                if msa_i in [5, 26, 11, 75, 48, 3, 61]:
                    to_sort_pos.append(def_list)
            # remove all breaks trends or trend no detected at all
            if data[2] == 1 or data[2] == -1:
                regardless_break.append(def_list)
                continue
            to_sort.append(def_list)

    samples = 35
    # sort by sequentiality, return 20 or less
    best_20 = sorted(to_sort, key=lambda x: (x[2], x[3]), reverse=True)[:samples]
    best_dict["Order"] = [i + 1 for i in range(len(best_20))]
    best_dict["Substitution"] = [x[0] for x in best_20]
    best_dict["Tree"] = [x[1] for x in best_20]
    best_dict["Sequentiality"] = [x[2] for x in best_20]
    best_dict["Fluctuation"] = [x[3] for x in best_20]
    best_dict["Break trend"] = [x[4] if x[4] != -1 else False for x in best_20]
    best_dict["N"] = [x[5] for x in best_20]
    best_dict["G"] = [x[6] for x in best_20]
    best_dict["AA CNT"] = [calc_aa_in_lineage(n) for n in best_dict["N"]]
    best_dict["Conservation Score"] = [conservations[x[7]] for x in best_20]

    best_20 = sorted(regardless_break, key=lambda x: (x[2]), reverse=True)[:samples]
    best_dict["Order"].extend([i + 1 for i in range(len(best_20))])
    best_dict["Substitution"].extend([x[0] for x in best_20])
    best_dict["Tree"].extend([x[1] for x in best_20])
    best_dict["Sequentiality"].extend([x[2] for x in best_20])
    best_dict["Fluctuation"].extend([x[3] for x in best_20])
    best_dict["Break trend"].extend([x[4] if x[4] != -1 else False for x in best_20])
    best_dict["N"].extend([x[5] for x in best_20])
    best_dict["G"].extend([x[6] for x in best_20])
    best_dict["AA CNT"].extend([calc_aa_in_lineage(n) for n in best_dict["N"][-samples:]])
    best_dict["Conservation Score"].extend([conservations[x[7]] for x in best_20])

    best_20 = sorted(regardless_break, key=lambda x: (x[3]), reverse=True)[:samples]
    best_dict["Order"].extend([i + 1 for i in range(len(best_20))])
    best_dict["Substitution"].extend([x[0] for x in best_20])
    best_dict["Tree"].extend([x[1] for x in best_20])
    best_dict["Sequentiality"].extend([x[2] for x in best_20])
    best_dict["Fluctuation"].extend([x[3] for x in best_20])
    best_dict["Break trend"].extend([x[4] if x[4] != -1 else False for x in best_20])
    best_dict["N"].extend([x[5] for x in best_20])
    best_dict["G"].extend([x[6] for x in best_20])
    best_dict["AA CNT"].extend([calc_aa_in_lineage(n) for n in best_dict["N"][-samples:]])
    best_dict["Conservation Score"].extend([conservations[x[7]] for x in best_20])

    if len(to_sort_pos) > 0:
        best_20 = sorted(to_sort_pos, key=lambda x: (x[2], x[3]), reverse=True)
        p0a8u6["Order_2"] = [i + 1 for i in range(len(best_20))]
        p0a8u6["Substitution_2"] = [x[0] for x in best_20]
        p0a8u6["Tree_2"] = [x[1] for x in best_20]
        p0a8u6["Sequentiality_2"] = [x[2] for x in best_20]
        p0a8u6["Fluctuation_2"] = [x[3] for x in best_20]
        p0a8u6["Break trend_2"] = [x[4] if x[4] != -1 else False for x in best_20]
        p0a8u6["N_2"] = [x[5] for x in best_20]
        p0a8u6["G_2"] = [x[6] for x in best_20]
        p0a8u6["AA CNT_2"] = [calc_aa_in_lineage(n) for n in p0a8u6["N_2"]]
        p0a8u6["Conservation Score"] = [conservations[x[7]] for x in best_20]
        # print(len(to_sort_pos))
    return best_dict, p0a8u6, indice_stats


def predict_level2(run: RunSetup):
    """Predicted level 2 consensus - now we have one sequence over all trees"""
    log_msg = f"Predicting Level 2 consensus for {run.protein_name} \n"

    input_pattern = os.path.join(run.tmp_files, 'tree[0-9][0-9][0-9][0-9]')
    trees = glob.glob(input_pattern)
    TREE_CNT = int(os.environ['tree_cnt'])

    data_path = os.environ['input']
    proteins = [run.protein_name]
    templ = "{}_set_{}/"

    # Fetch all necessary data
    with open(os.path.join(run.index_fld, "aa_idx.pkl"), "rb") as f:
        aa_to_idx = pickle.load(f)

    with open(os.path.join(run.index_fld, "file_names.pkl"), "rb") as f:
        file_names = pickle.load(f)

    with open(os.path.join(run.index_fld, "aa_indice_names.pkl"), "rb") as f:
        selected_indices = pickle.load(f)

    aa_ind_matrixes = dict()
    aa_idx_names = list()
    for file in file_names:
        aa_idx_names.append(file.split(".")[0])
        with open(os.path.join(run.index_fld, file), "rb") as f:
            aa_ind_matrixes[selected_indices[len(aa_idx_names) - 1]] = pickle.load(f)

    CONFIDENCE_LEVEL = run.confidence_level
    # prepare csv files for every
    for protein in proteins:
        gen_dir = f"outputs/{protein}/level2/" os.path.join()
        os.makedirs(gen_dir, exist_ok=True)
        tmp_dir = os.path.join(run.tmp_files, "level2")  # f"outputs/{protein}/level2/script_tmp/"
        os.makedirs(tmp_dir, exist_ok=True)

        log_msg += "CONFIDENCE LEVEL {}".format(CONFIDENCE_LEVEL)+"\n"
        print(log_msg)

        # fetch all r2 and slope metrics
        metrics = {"r2": {f: [] for f in selected_indices},
                   "slope": {f: [] for f in selected_indices},
                   "sequence_trend": {f: [] for f in selected_indices},
                   "sequence_trend_break": {f: [] for f in selected_indices},
                   "fluctuations": {f: [] for f in selected_indices},
                   "values": {f: [] for f in selected_indices}
                   }

        indicis_sequences = dict()
        for ind_name in selected_indices:
            indicis_sequences[ind_name] = list()

        # catch all sequences for protein, and all statistics from tmp files
        for folder_i, tree_fld in enumerate(trees):
            indices_file = open(os.path.join(tree_fld, "prediction_indices_{}.fasta".format(CONFIDENCE_LEVEL)))
            # TODO check indices sequences what does it mean
            id_name = ""
            for line in indices_file.readlines():
                line = line.strip()
                if line == "" or line == "\n":
                    continue
                if line[0] == ">":
                    id_name = line[1:]
                    continue
                indicis_sequences[id_name].append(line)
            indices_file.close()

            for ind in selected_indices:
                with open(os.path.join(tree_fld, f"r2_{ind}.pkl"), "rb") as f:
                    metrics["r2"][ind].append(pickle.load(f))
                with open(os.path.join(tree_fld, f"slope_{ind}.pkl"), "rb") as f:
                    metrics["slope"][ind].append(pickle.load(f))
                with open(os.path.join(tree_fld, f"sequentiality_{ind}.pkl"), "rb") as f:
                    metrics["sequence_trend"][ind].append(pickle.load(f))
                with open(os.path.join(tree_fld, f"sequentiality2_{ind}.pkl"), "rb") as f:
                    metrics["sequence_trend_break"][ind].append(pickle.load(f))
                with open(os.path.join(tree_fld, f"fluctuation_{ind}.pkl"), "rb") as f:
                    metrics["fluctuations"][ind].append(pickle.load(f))
                with open(os.path.join(tree_fld, f"values_{ind}.pkl"), "rb") as f:
                    metrics["values"][ind].append(pickle.load(f))

        # Prepare WT
        wt_dict, _ = load_msa(str(os.path.join(run.ground_truth, run.protein_name)))
        wt = wt_dict[run.query]

        indices_stats = {}
        # create aa indices fasta for alignment and align them
        for aa_id in selected_indices:
            print(aa_id)
            indices_stats[aa_id] = []
            wt_fasta_name = aa_id + "_{}_toAlign.fasta".format(CONFIDENCE_LEVEL)
            of = open(tmp_dir + wt_fasta_name, "w")
            for seq_i, seq in enumerate(indicis_sequences[aa_id]):
                of.write(">seq_tree_{}\n".format(seq_i))
                of.write(seq.replace("-", ""))
                of.write("\n")
            of.close()

            msa_file = "{}_{}_msa.fasta".format(aa_id, CONFIDENCE_LEVEL)
            os.system("clustalo -i " + tmp_dir + wt_fasta_name + " -o " + tmp_dir + msa_file + " --force")

            # create statistics for indices
            mf = open(tmp_dir + msa_file)
            alignment = dict()
            pred_cnt, loading_sequence = 0, ""
            for line in mf.readlines():
                line = line.strip()
                if line == "" or line == "\n":
                    continue
                if line[0] == ">":
                    if pred_cnt == 1:
                        alignment[seq_name] = loading_sequence
                        loading_sequence = ""
                    seq_name = line[1:]
                    pred_cnt = 1
                    continue
                loading_sequence += line
            mf.close()
            # Last sequence not in
            alignment[seq_name] = loading_sequence

            wt_fasta_file_path = os.path.join(run.fasta, f"{aa_id}_{CONFIDENCE_LEVEL}.fasta")
            msa_file_path = os.path.join(run.fasta, f"{aa_id}_{CONFIDENCE_LEVEL}_msa.fasta")
            alignment = clustalo_sequences(wt_fasta_file_path, msa_file_path, )

            # create mapping of MSA position to original one
            seq_mapping = dict()
            trees_map_to_msa = list()
            over_tree_msa_len = len(list(alignment.values())[0])
            for name, seq in alignment.items():
                seq_mapping[name] = [i for i, aa in enumerate(seq) if aa != '-']  # items are indexes to msa
                trees_map_to_msa.append(seq_mapping[name])

            align_len = len(list(alignment.items())[0][1])
            stats = dict()
            for pos in range(align_len):
                stats[str(pos)] = {aa: 0 for aa in aa_to_idx}
                stats[str(pos)]["-"] = 0
            # Count occurrences
            for seq in alignment.values():
                for pos_i, pos in enumerate(seq):
                    stats[str(pos_i)][pos] += 1

            result_dict = dict()
            result_dict["total"] = TREE_CNT
            predicted_sequence = ""
            for pos, vals in stats.items():
                ordered_level = {k: v for k, v in sorted(vals.items(), key=lambda item: item[1])}
                result_dict[pos] = list(list(ordered_level.items())[-1])
                predicted_sequence += list(ordered_level.items())[-1][0]
            df = pd.DataFrame.from_dict(result_dict)
            df.to_csv(tmp_dir + "{}_average_{}.csv".format(aa_id, CONFIDENCE_LEVEL))

            # compare with WT
            wt_fasta_name = aa_id + "_{}_toAlignWithWT.fasta".format(CONFIDENCE_LEVEL)
            of = open(tmp_dir + wt_fasta_name, "w")

            of.write(">wt_{}\n".format(protein))
            of.write(wt.replace("-", ""))
            of.write(">predicted_{}\n".format(aa_id))
            of.write(predicted_sequence.replace("-", ""))
            of.write("\n")
            of.close()

            msa_comparison = "{}_{}_comparison.fasta".format(aa_id, CONFIDENCE_LEVEL)
            os.system("clustalo -i " + tmp_dir + wt_fasta_name + " -o " + tmp_dir + msa_comparison + " --force")

            # create statistics for indices
            mf = open(tmp_dir + msa_comparison)
            alignment_comp = dict()
            loading_sequence = ""
            for line in mf.readlines():
                line = line.strip()
                if line == "" or line == "\n":
                    if loading_sequence != "":
                        alignment_comp[seq_name] = loading_sequence
                    continue
                if line[0] == ">":
                    if loading_sequence == "":
                        seq_name = line[1:]
                        continue
                    alignment_comp[seq_name] = loading_sequence
                    loading_sequence = ""
                    seq_name = line[1:]
                    continue
                loading_sequence += line
            mf.close()
            # Last sequence not in
            alignment_comp[seq_name] = loading_sequence

            comparison_dict = {"Notes": [
                "{}".format(protein),
                "prediction",
                "{} index similarity".format(aa_id)]}

            wt_al, pred = list(alignment_comp.values())  # just 2 sequence in there
            for pos_i, (w, p) in enumerate(zip(wt_al, pred)):
                if w == "-" or p == "-":
                    d = ""
                else:
                    d = aa_ind_matrixes[aa_id][aa_to_idx[w]][aa_to_idx[p]]
                comparison_dict[pos_i] = [w, p, d]

            # metric statistics
            metrics_dict = {"Notes": ['Sequence',
                                      f'{protein}', 'consensus freq',
                                      'mean r2', 'std r2', 'r2 max min', 'r2 N',
                                      'mean slope', 'std slope', 'slope max min',
                                      'slope N', ""]}

            metrics_dict["Notes"].extend(['tree{}, (R,Sl,AA,S,F,B,P,G)'.format(t_n) for t_n in range(TREE_CNT)])

            metrics_stats_r2, r2_col_vals, voting_predictions = create_stats_with_metric2(metrics["r2"][aa_id],
                                                                                          alignment, trees_map_to_msa,
                                                                                          over_tree_msa_len)
            # voting is same for
            metrics_stats_slope, col_vals, _ = create_stats_with_metric2(metrics["slope"][aa_id], alignment,
                                                                         trees_map_to_msa, over_tree_msa_len)

            trends = sequence_trend_metrics(metrics, aa_id, trees_map_to_msa, over_tree_msa_len)
            aligned_wt = ""
            no_gaps_idx = 0
            for pos_i in range(over_tree_msa_len):
                seq_aa = predicted_sequence[pos_i]
                wt_res = ''
                if seq_aa != '-':
                    wt_res = wt_al[no_gaps_idx]
                    no_gaps_idx += 1

                aligned_wt += "-" if wt_res == '' else wt_res

                voted_aa, voted_freq = voting_by_prediction(voting_predictions[pos_i])

                col_trends = trends[pos_i]
                mean, met_stddev, met_max, met_min, cnt = metrics_stats_r2[pos_i]
                mean_s, met_stddev_s, met_max_s, met_min_s, cnt_s = metrics_stats_slope[pos_i]
                col_csv_val = [seq_aa, wt_res, f"{voted_aa}, {voted_freq}",
                               "{:.3f}".format(mean), "{:.3f}".format(met_stddev),
                               "{:.3f}, {:.3f}".format(met_max, met_min), cnt,
                               "{:.3f}".format(mean_s), "{:.3f}".format(met_stddev_s),
                               "{:.3f}, {:.3f}".format(met_max_s, met_min_s), cnt_s,
                               ""]

                cons_score = None

                tree_records = ["{:.3f}, {:.3f}, {}".format(r2[0], sl[0], sl[1]) for r2, sl in zip(r2_col_vals[pos_i],
                                                                                                   col_vals[pos_i])]
                for trend_i in range(TREE_CNT):
                    trend = col_trends[trend_i]
                    tree_records[trend_i] += f", {trend[0]}, {trend[1]}, {trend[2]}, {trend[3]}, {trend[4]}, " \
                                             f"plot_col_{trend[5]}"
                col_csv_val.extend(tree_records)
                metrics_dict[pos_i] = col_csv_val

            best_met, positions, aa_stats = detect_best_metric(trends, alignment, aligned_wt, len(metrics_dict[0]), protein,
                                                               cons_score)
            indices_stats[aa_id] = aa_stats

            df = pd.DataFrame.from_dict(comparison_dict)
            df.to_csv(tmp_dir + "{}_averageWT_{}.csv".format(aa_id, CONFIDENCE_LEVEL))

            # append the best metrics dict to csv
            metrics_dict.update(best_met)

            if protein == "P0A8U6":
                metrics_dict.update(positions)

            max_len = max([len(x) for _, x in metrics_dict.items()])

            # add lines to dictionary to solve error of different line count
            for k, ls in metrics_dict.items():
                rows_missing = max_len - len(ls)
                filling = ["" for _ in range(rows_missing)]
                ls.extend(filling)

            df = pd.DataFrame.from_dict(metrics_dict)
            df.to_csv(gen_dir + "{}_metrics_{}.csv".format(aa_id, CONFIDENCE_LEVEL))

            # store final
            final_seq = "{}_FINAL_SEQUENCE_{}.fasta".format(aa_id, CONFIDENCE_LEVEL)
            mf = open(gen_dir + final_seq, "w")
            mf.write(pred.replace("-", ""))
            mf.close()

        # fill up stats with blank cells
        max_len = max([len(x) for _, x in indices_stats.items()])
        for k, s in indices_stats.items():
            rows_missing = max_len - len(s)
            filling = ["" for _ in range(rows_missing)]
            s.extend(filling)

        df = pd.DataFrame.from_dict(indices_stats)
        df.to_csv(tmp_dir + "indices_transitions{}.csv".format(CONFIDENCE_LEVEL))
        with open(tmp_dir + "indices_transitions{}.pkl".format(CONFIDENCE_LEVEL), "wb") as f:
            pickle.dump(indices_stats, f)
