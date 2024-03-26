# SSP - Successor Sequence Predictor

This repository of scripts aims to investigate the evolutionary successor reconstruction strategy (SSP)
as a complementary method to the Ancestral Sequence reconstruction (ASR). We evaluate SSP for a design of enhanced 
variance in the evolutionary-based scheme according to physic-chemical AA properties of multiple AA indices. 

The scripts look for linear trends in evolution (simulated by ancestral sequences in a phylogenetic tree) 
over a selected set of AA indices to get multiple suggested AA-index related mutations. 

We conducted an analysis on a specific protein using multiple randomly generated phylogenetic trees,
including ancestors.
Firstly, we made a prediction for each tree and amino acid index (level 1).
Secondly, we established a consensus at level 2 based on the predictions from all the phylogenetic trees.
Finally, we suggested a level 3 sequence by combining all the level 2 predictions.

## Installation

Please install clustalo to your system simply running `sudo apt-get -y install clustalo` and check it is accessible in terminal 
as `clustalo -h`.

**These instructions assume you have the python package manager `conda` installed.**

Go to the directory where you wish to have your project be stored and run
```
cd folder/you/wish/to/place/your/project
git clone git@git.loschmidt.cz/ssp
cd ssp
```
Create conda environment for the project
```
conda create --name ssp python=3.10
```
Activate environment and install requirement and make our project visible in the environment
```
conda activate ssp
conda install --file requirements.txt
pip install -e .  # by this you make our successor library visible in the environment
```
## Configuration
`Yaml` configuration file is required when running experiments. To make setup of experiments more 
cloud friendly the `cli/config.py` client and `configs/config_template.yaml` exist.

Copy template file:
```
cp successors/configs/config_template.yaml my_cofig.yaml
```
Set all fields by current experiment requirements, e.g. set out dir
All supported options and their descriptions for a running pipeline can be found there.

It is worth mentioning, that `transition` option set to **YES** can case 'no mutational design' 
for shallow phylogenetic trees (algorithm then applies prediction for evolutionary position with three or more 
transitions as we consider this as a strong evolutionary signal). 

For deep trees, `transition` option set to **NO** can cause to noisy sequence generation for deep trees (seven and more node from root to leaf).

#### Dataset description

Following supplement S2 manual for creating dataset for SSP prediction you can create data for your analysis
using SigClust and Blast tools. 
To replicate **step 3 from Supplement 2** and prepare data batches for phylogenetic tree construction using 
FireProtAsr web tool use script `create_sets.py`. You can run the script as following:
```
python successors/create_sets.py --input_seq my_blast_sequences.fasta --clusters SigClust_ouput --query my_query_id --output_dir my_SPP_analysis 
```

Scripts account for the dataset folder structure, set in `input` config parameter, as following:
1. each tree has a separate folder (tree1, tree2 ...) 
2. treeN folder must have **msa.fasta** and **ancestralTree.tree** files. (keep names as these)
   1. **msa.fasta**: whole tree MSA or just sequences on the tree trajectory from the root to the 'query'
   2. **ancestralTree.tree**: newick format of the phylogenetic tree

#### Adding custom AA Index
To add custom index to the study, one must append AA-wise index values to `successors/indices/aaindex.csv` and `successors/aa_indexes.py` files.

## Running

To run the pipeline, you will find results in the corresponding folder specified in config file

For demonstration,
you can run example with LGK protein and example config file in `results/lgk/config/lgk_config.yaml` file.

### Prediction
Make level1 prediction (successors per tree and AA-index):
```
python path/to/the/main/cli/main.py predict /path/to/the/conf.yaml
```
All processing information are logged in the experiment directory folder  `experiment_path_dir/logs/log_file.txt`

### AA-INDEX sequence consensus
To get level2 prediction (successor per AA-INDEX)
```
python path/to/the/main/cli/main.py level2 /path/to/the/conf.yaml
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
python path/to/the/main/cli/main.py level3 /path/to/the/conf.yaml 
```
These results might be of particular interest as there is a complex overview over all mutations and their AA-indices, 
so one can examine why the given suggestion appeared.

##### Final consensus sequence
The final consensus over all indices and trees can be found in `/path/to/output/dir/results/final_majority_consensus.fasta` 
along with a wild type sequence for comparison. 
Please note
that a final consensus sequence may vary over multiple iterations
in case there are two similar AA frequencies for particular position. 



