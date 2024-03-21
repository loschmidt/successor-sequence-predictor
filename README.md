# SSP - Successor Sequence Predictor

This Python package aims to investigate the evolutionary successor reconstruction strategy (SSP)
as a complementary method to the Ancestral Sequence reconstruction (ASR). We evaluate SSP for a design of enhanced 
variance in the evolutionary-based scheme according to physic-chemical amino-acid (AA) properties of multiple AA indices. 

The package looks for linear trends in evolution (simulated by ancestral sequences in a phylogenetic tree) 
over a selected set of AA indices to get multiple suggested AA-index related mutations. 

We conducted an analysis on a specific protein using multiple randomly generated phylogenetic trees,
including ancestors.
Firstly, we made a prediction for each tree and amino acid index (level 1).
Secondly, we established a consensus at level 2 based on the predictions from all the phylogenetic trees.
Finally, we suggested a level 3 sequence by combining all the level 2 predictions.

## Prerequisites

- Clustal Omega
- Python >= 3.10

If you use Debian or Ubuntu based Linux distribution, you can install Clustal Omega from software repository:
```bash
sudo apt-get install clustalo
```

For other systems, follow instructions on the Clustal Omega download page: http://www.clustal.org/omega/

To verify that Clustal Omega is installed, run: `clustalo -h`.

## Installation

Install SSP package using pip:
```
pip3 install https://github.com/loschmidt/successor-sequence-predictor/archive/refs/heads/main.zip
```

Verify that SSP package is installed:
```
loschmidt-ssp --help
```

## Configuration
YAML configuration file is required when running experiments. To make setup of experiments more cloud friendly, the `cli/config.py` client and `configs/config_template.yaml` exist.

Copy template file:
```bash
cp successors/configs/config_template.yaml my_cofig.yaml
```
Set all fields by current experiment requirements, e.g. set out dir
All supported options and their descriptions for a running pipeline can be found there.

It is worth mentioning, that `transition` option set to **YES** can case 'no mutational design' 
for shallow phylogenetic trees (algorithm then applies prediction for evolutionary position with three or more 
transitions as we consider this as a strong evolutionary signal). 

For deep trees, `transition` option set to **NO** can cause to noisy sequence generation for deep trees (seven and more node from root to leaf).

#### Dataset description
Scripts account for the dataset folder structure, set in `input` config parameter, as following:
1. each tree has a separate folder (tree1, tree2 ...) 
2. treeN folder must have **msa.fasta** and **ancestralTree.tree** files. (keep names as these)
   1. **msa.fasta**: whole tree MSA or just sequences on the tree trajectory from the root to the 'query'
   2. **ancestralTree.tree**: newick format of the phylogenetic tree

#### Adding custom AA Index
To add custom index to the study, one must append AA-wise index values to `successors/indices/aaindex.csv` and `successors/aa_indexes.py` files.

## Running

To run the pipeline, you will find results in the corresponding folder specified in the config file.

For demonstration, you can run example with LGK protein and example config file in `results/lgk/config/lgk_config.yaml` file.

### Prediction
Make level1 prediction (successors per tree and AA-index):
```
loschmidt-ssp predict /path/to/the/conf.yaml
```
All processing information are logged in the experiment directory folder  `experiment_path_dir/logs/log_file.txt`

### AA-INDEX sequence consensus
To get level2 prediction (successor per AA-INDEX)
```
loschmidt-ssp level2 /path/to/the/conf.yaml
```
This command will generate quite many files in `results/level2` folder which you can investigate, 
but it is not necessary. It is worth mentioning that this folder includes csv files *AA-INDEX*_averageWT_0.csv with 
AA-index similarities of replaced AAs in the final index consensus (1 - no change in property)

### Final Consensus
Generate level3 statistics for individual substitutions in the tree and indices. 
It generates two statistical files:
1. **a3_*protein*_metrics_0.csv** - list and statistics for exact match of substitutions (same position and replaced AA)
2. **b3_*protein*_metrics_0.csv** - list and statistics for substitutions from the same position but can be replaced by any AA  
```
loschmidt-ssp level3 /path/to/the/conf.yaml 
```
These results might be of particular interest as there is a complex overview over all mutations and their AA-indices, 
so one can examine why the given suggestion appeared.

##### Final consensus sequence
The final consensus over all indices and trees can be found in `/path/to/output/dir/results/final_majority_consensus.fasta`  along with a wild type sequence for comparison. 
Please note that a final consensus sequence may vary over multiple iterations
in case there are two similar AA frequencies for particular position. 
