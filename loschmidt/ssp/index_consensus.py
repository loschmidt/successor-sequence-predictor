__author__ = "Pavel Kohout xkohou15@vutbr.cz"
__date__ = "2024/02/19"
__description__ = "For every AA indice make an average and compute WT similarity over all trees"

import glob
import math
import os
import pickle
import statistics

import pandas as pd

from loschmidt.ssp.parser_handler import RunSetup
from loschmidt.ssp.utils import load_msa, clustalo_sequences


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
    """
    Compute the metrics in the evolutionary trajectory
    :param metrics:
    :param feature:
    :param tree_map_to_msa:
    :param indice_msa_len:
    :return:
    """
    seq_data = metrics["sequence_trend"][feature]
    break_trend = metrics["sequence_trend_break"][feature]
    fluctuation = metrics["fluctuations"][feature]
    ys = metrics["values"][feature]

    templates = []
    for msa_pos in range(indice_msa_len):
        # S,F,B,P,G, plot
        template = [[0, 0, -1, 0, 0, -1] for _ in range(len(tree_map_to_msa))]
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


def detect_best_metric(trends, trees_msa, aligned_wt, dict_line_cnt, protein, run):
    def calc_aa_in_lineage(n):
        aa_cnt, cumsum = 0, 0
        while True:
            cumsum += aa_cnt
            if cumsum == n:
                return aa_cnt
            aa_cnt += 1

    # conservation
    conservations = [0 for _ in range(1000)]  # fill it to make code clear, allocate 1000 columns in csv
    # for visualization purposes if alignment is not well aligned at the beginning just for visualization no for prediction
    if run.conservation_offset:
        conservations = ['' for _ in range(int(run.conservation_offset))] + conservations

    # convert MSA dict to array
    msa_index_list = list()
    for _, sequence in trees_msa.items():
        msa_index_list.append(sequence)

    best_dict = {"Individul trees": ["" for _ in range(dict_line_cnt)],
                 "Order": [], "Substitution": [], "Tree": [], "Sequentiality": [], "Fluctuation": [], "Break trend": [],
                 "N": [], "G": [], "AA CNT": [], "Conservation Score": []}
    custom_selection = {f"Custom {protein}": ["" for _ in range(dict_line_cnt)],
              "Order_2": [], "Substitution_2": [], "Tree_2": [], "Sequentiality_2": [], "Fluctuation_2": [],
              "Break trend_2": [],
              "N_2": [], "G_2": [], "AA CNT_2": [], "Conservation Score": []}
    # S,F,B,P,G
    to_sort = list()
    regardless_break = list()
    indice_stats = list()
    custom_selection_positions = list()
    for msa_i, trees_data in enumerate(trends):
        for tree_i, data in enumerate(trees_data):
            # substitutions, tree
            if aligned_wt[msa_i] == msa_index_list[tree_i][msa_i]:
                continue  # that is not a substitution
            def_list = [f"{aligned_wt[msa_i]}{msa_i}{msa_index_list[tree_i][msa_i]}", tree_i]
            def_list.extend(data)
            def_list.extend([msa_i])
            indice_stats.append(def_list)
            # positions targets for custom inspection
            if msa_i in run.highlight_pos:  #[5, 26, 11, 75, 48, 3, 61]
                custom_selection_positions.append(def_list)
            # remove all breaks trends or trend no detected at all
            if data[2] == 1 or data[2] == -1:
                regardless_break.append(def_list)
                continue
            to_sort.append(def_list)

    samples = 20
    # sort by sequentiality, return 20 or less
    best_n = sorted(to_sort, key=lambda x: (x[2], x[3]), reverse=True)[:samples]
    best_dict["Order"] = [i + 1 for i in range(len(best_n))]
    best_dict["Substitution"] = [x[0] for x in best_n]
    best_dict["Tree"] = [x[1] for x in best_n]
    best_dict["Sequentiality"] = [x[2] for x in best_n]
    best_dict["Fluctuation"] = [x[3] for x in best_n]
    best_dict["Break trend"] = [x[4] if x[4] != -1 else False for x in best_n]
    best_dict["N"] = [x[5] for x in best_n]
    best_dict["G"] = [x[6] for x in best_n]
    best_dict["AA CNT"] = [calc_aa_in_lineage(n) for n in best_dict["N"]]
    best_dict["Conservation Score"] = [conservations[x[7]] for x in best_n]

    best_n = sorted(regardless_break, key=lambda x: (x[2]), reverse=True)[:samples]
    best_dict["Order"].extend([i + 1 for i in range(len(best_n))])
    best_dict["Substitution"].extend([x[0] for x in best_n])
    best_dict["Tree"].extend([x[1] for x in best_n])
    best_dict["Sequentiality"].extend([x[2] for x in best_n])
    best_dict["Fluctuation"].extend([x[3] for x in best_n])
    best_dict["Break trend"].extend([x[4] if x[4] != -1 else False for x in best_n])
    best_dict["N"].extend([x[5] for x in best_n])
    best_dict["G"].extend([x[6] for x in best_n])
    best_dict["AA CNT"].extend([calc_aa_in_lineage(n) for n in best_dict["N"][-samples:]])
    best_dict["Conservation Score"].extend([conservations[x[7]] for x in best_n])

    best_n = sorted(regardless_break, key=lambda x: (x[3]), reverse=True)[:samples]
    best_dict["Order"].extend([i + 1 for i in range(len(best_n))])
    best_dict["Substitution"].extend([x[0] for x in best_n])
    best_dict["Tree"].extend([x[1] for x in best_n])
    best_dict["Sequentiality"].extend([x[2] for x in best_n])
    best_dict["Fluctuation"].extend([x[3] for x in best_n])
    best_dict["Break trend"].extend([x[4] if x[4] != -1 else False for x in best_n])
    best_dict["N"].extend([x[5] for x in best_n])
    best_dict["G"].extend([x[6] for x in best_n])
    best_dict["AA CNT"].extend([calc_aa_in_lineage(n) for n in best_dict["N"][-samples:]])
    best_dict["Conservation Score"].extend([conservations[x[7]] for x in best_n])

    # THIS IS CUSTOM HIGHLIGHT FOR SPECIAL POSITIONS WE ARE INTERESTED IN
    if len(custom_selection_positions) > 0:
        best_n = sorted(custom_selection_positions, key=lambda x: (x[2], x[3]), reverse=True)
        custom_selection["Order_2"] = [i + 1 for i in range(len(best_n))]
        custom_selection["Substitution_2"] = [x[0] for x in best_n]
        custom_selection["Tree_2"] = [x[1] for x in best_n]
        custom_selection["Sequentiality_2"] = [x[2] for x in best_n]
        custom_selection["Fluctuation_2"] = [x[3] for x in best_n]
        custom_selection["Break trend_2"] = [x[4] if x[4] != -1 else False for x in best_n]
        custom_selection["N_2"] = [x[5] for x in best_n]
        custom_selection["G_2"] = [x[6] for x in best_n]
        custom_selection["AA CNT_2"] = [calc_aa_in_lineage(n) for n in custom_selection["N_2"]]
        custom_selection["Conservation Score"] = [conservations[x[7]] for x in best_n]
        # print(len(to_sort_pos))
    return best_dict, custom_selection, indice_stats


def predict_level2(run: RunSetup):
    """Predicted level 2 consensus - now we have one sequence over all trees"""
    log_msg = f"Predicting Level 2 consensus for {run.protein_name} \n"

    input_pattern = os.path.join(run.tmp_files, 'tree*')
    trees = glob.glob(input_pattern)
    tree_cnt = len(trees)  # int(os.environ['tree_cnt'])
    protein = [run.protein_name]
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
    level2_dir = os.path.join(run.results, "level2")
    os.makedirs(level2_dir, exist_ok=True)

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
    for tree_fld in trees:
        # load predicted index sequences for all trees and sort them to one variable by index
        indices_file_path = os.path.join(tree_fld, "prediction_indices_{}.fasta".format(CONFIDENCE_LEVEL))
        indices_predictions_per_tree, _ = load_msa(indices_file_path)
        for aa_idx, seq in indices_predictions_per_tree.items():
            indicis_sequences[aa_idx].append(seq)
        # get calculate statistics over all trees
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
        print("Processing... " + aa_id)
        indices_stats[aa_id] = []

        # align all predictions for one AA index over all trees
        wt_fasta_file_path = os.path.join(run.fasta, f"{aa_id}_{CONFIDENCE_LEVEL}.fasta")
        msa_file_path = os.path.join(run.fasta, f"{aa_id}_{CONFIDENCE_LEVEL}_msa.fasta")
        sequences_to_align = {f"seq_tree_{i_tree}": s for i_tree, s in enumerate(indicis_sequences[aa_id])}
        alignment, msa_len = clustalo_sequences(wt_fasta_file_path, msa_file_path, sequences_to_align)

        # create mapping of MSA position to original one
        seq_mapping = dict()
        trees_map_to_msa = list()
        for name, seq in alignment.items():
            seq_mapping[name] = [i for i, aa in enumerate(seq) if aa != '-']  # items are indexes to msa
            trees_map_to_msa.append(seq_mapping[name])

        per_position_aa_frequencies = dict()
        for pos in range(msa_len):
            per_position_aa_frequencies[str(pos)] = {aa: 0 for aa in aa_to_idx}
            per_position_aa_frequencies[str(pos)]["-"] = 0
        # Count occurrences
        for seq in alignment.values():
            for pos_i, pos in enumerate(seq):
                per_position_aa_frequencies[str(pos_i)][pos] += 1

        result_dict = dict()
        result_dict["total"] = tree_cnt
        predicted_sequence = ""
        for pos, vals in per_position_aa_frequencies.items():
            ordered_level = {k: v for k, v in sorted(vals.items(), key=lambda item: item[1])}
            result_dict[pos] = list(list(ordered_level.items())[-1])
            predicted_sequence += list(ordered_level.items())[-1][0]
        df = pd.DataFrame.from_dict(result_dict)
        df.to_csv(os.path.join(level2_dir, f"{aa_id}_average_{CONFIDENCE_LEVEL}.csv"))

        # compare a predicted AA index sequence with original WT
        wt_predicted_fasta_path = os.path.join(run.fasta, f"{aa_id}_{CONFIDENCE_LEVEL}_toAlignWithWT.fasta")
        msa_wt_comparison = os.path.join(run.fasta, "{}_{}_comparison.fasta".format(aa_id, CONFIDENCE_LEVEL))

        wt_predicted_fasta = {
            f"wt_{protein}": wt.replace("_", ""),
            f"predicted_{aa_id}": predicted_sequence.replace("-", "")
        }

        alignment_comp, comp_msa_len = clustalo_sequences(wt_predicted_fasta_path, msa_wt_comparison, wt_predicted_fasta)

        ######################################################
        # Now proceed logic for generation of a statistic report with sequentiality fluctuation per AA aligned to wt
        #   highly customized report, please modify to your needs

        comparison_dict = {"Notes": [
            "{}".format(protein),
            "prediction",
            "{} index similarity".format(aa_id)]}

        # get similarity of given index wt AA to successor AA
        wt_al, pred = list(alignment_comp.values())  # just 2 sequence in there
        for pos_i, (w, p) in enumerate(zip(wt_al, pred)):
            if w == "-" or p == "-":
                similarity = ""
            else:
                similarity = aa_ind_matrixes[aa_id][aa_to_idx[w]][aa_to_idx[p]]
            comparison_dict[pos_i] = [w, p, similarity]

        # metric statistics
        metrics_dict = {"Notes": ['Sequence',
                                  f'{protein}', 'consensus freq',
                                  'mean r2', 'std r2', 'r2 max min', 'r2 N',
                                  'mean slope', 'std slope', 'slope max min',
                                  'slope N', ""]}

        metrics_dict["Notes"].extend(['tree{}, (R,Sl,AA,S,F,B,P,G)'.format(t_n) for t_n in range(tree_cnt)])

        metrics_stats_r2, r2_col_vals, voting_predictions = create_stats_with_metric2(metrics["r2"][aa_id],
                                                                                      alignment, trees_map_to_msa,
                                                                                      msa_len)
        # voting is the same for
        metrics_stats_slope, col_vals, _ = create_stats_with_metric2(metrics["slope"][aa_id], alignment,
                                                                     trees_map_to_msa, msa_len)

        trends = sequence_trend_metrics(metrics, aa_id, trees_map_to_msa, msa_len)

        # Create the final csv file per AA index with aligned maximum voted prediction and all statistics together
        aligned_wt = ""
        no_gaps_idx = 0
        for pos_i in range(msa_len):
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

            tree_records = ["{:.3f}, {:.3f}, {}".format(r2[0], sl[0], sl[1]) for r2, sl in zip(r2_col_vals[pos_i],
                                                                                               col_vals[pos_i])]
            for trend_i in range(tree_cnt):
                trend = col_trends[trend_i]
                tree_records[trend_i] += f", {trend[0]}, {trend[1]}, {trend[2]}, {trend[3]}, {trend[4]}, " \
                                         f"plot_col_{trend[5]}"
            col_csv_val.extend(tree_records)
            metrics_dict[pos_i] = col_csv_val

        best_met, positions, aa_stats = detect_best_metric(trends, alignment, aligned_wt, len(metrics_dict[0]), protein,
                                                           run)
        indices_stats[aa_id] = aa_stats

        df = pd.DataFrame.from_dict(comparison_dict)
        df.to_csv(os.path.join(level2_dir, "{}_averageWT_{}.csv".format(aa_id, CONFIDENCE_LEVEL)))

        # append the best metrics dict to csv
        metrics_dict.update(best_met)

        # in the case of used custom highlight show also that one
        if len(run.highlight_pos):
            metrics_dict.update(positions)

        max_len = max([len(x) for _, x in metrics_dict.items()])

        # add lines to dictionary to solve error of different line count
        for k, ls in metrics_dict.items():
            rows_missing = max_len - len(ls)
            filling = ["" for _ in range(rows_missing)]
            ls.extend(filling)

        df = pd.DataFrame.from_dict(metrics_dict)
        report_metric_path = os.path.join(level2_dir, f"{aa_id}_metrics_{CONFIDENCE_LEVEL}.csv")
        # df.to_csv(gen_dir + "{}_metrics_{}.csv".format(aa_id, CONFIDENCE_LEVEL))
        df.to_csv(report_metric_path)

        # store final
        # final_seq = "{}_FINAL_SEQUENCE_{}.fasta".format(aa_id, CONFIDENCE_LEVEL)
        predicated_consensus_path = os.path.join(level2_dir, f"{aa_id}_FINAL_SEQUENCE_{CONFIDENCE_LEVEL}.fasta")
        mf = open(predicated_consensus_path, "w")
        mf.write(pred.replace("-", ""))
        mf.close()

    # fill up stats with blank cells
    max_len = max([len(x) for _, x in indices_stats.items()])
    for k, s in indices_stats.items():
        rows_missing = max_len - len(s)
        filling = ["" for _ in range(rows_missing)]
        s.extend(filling)

    df = pd.DataFrame.from_dict(indices_stats)
    indices_stats_csv_path = os.path.join(level2_dir, f"indices_transitions{CONFIDENCE_LEVEL}.csv")
    indices_stats_pkl_path = os.path.join(level2_dir, f"indices_transitions{CONFIDENCE_LEVEL}.pkl")
    df.to_csv(indices_stats_csv_path)
    with open(indices_stats_pkl_path, "wb") as f:
        pickle.dump(indices_stats, f)
