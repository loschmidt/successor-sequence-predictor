import os
import pickle

from successors.parser_handler import RunSetup
from successors.utils import load_msa


def majority_voter(run: RunSetup):
    protein = run.protein_name

    # Fetch all necessary data
    with open(os.path.join(run.index_fld, "aa_indice_names.pkl"), "rb") as f:
        selected_indices = pickle.load(f)

    CONFIDENCE_LEVEL = run.confidence_level

    wt_seq = ""
    all_predictions = []
    final_prediction = ""
    max_len = 0
    for i, aa_index in enumerate(selected_indices):
        comparison_file = os.path.join(run.fasta, f"{aa_index}_{CONFIDENCE_LEVEL}_comparison.fasta")
        fasta_seq, alignment_len = load_msa(comparison_file)
        if i == 0:
            wt_seq = fasta_seq[list(fasta_seq.keys())[0]]
        all_predictions.append(fasta_seq[list(fasta_seq.keys())[1]])
        if alignment_len > max_len:
            max_len = alignment_len

    for i in range(max_len):
        pos_pred = [pred[i] for pred in all_predictions]
        final_aa = max(set(pos_pred), key=pos_pred.count)
        final_prediction += final_aa

    final_consensus_path = os.path.join(run.results, "final_majority_consensus.fasta")
    with open(final_consensus_path, "w") as f:
        for k, v in {f">wt_{protein}": wt_seq, ">final_pred": final_prediction}.items():
            f.write(f"{k}\n")
            for ch in [v[i:i+60] for i in range(0, len(v), 60)]:
                f.write("".join(ch) + "\n")
            f.write("\n")
    print(f"  Consensus was written in {final_consensus_path}")
