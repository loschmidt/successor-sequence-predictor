# Vaencestors

This project attempts to serve as backend for generation of ancestral like sequences

Scripts in this repository are expected to be used as the backend service for 
FireProt tool

## Installation

Please install clustalo to your system simply running `sudo apt-get -y install clustalo`

**These instructions assume you have the python package manager `conda` installed.**

Go to the directory where you wish to have your project be stored and run
```
cd folder/you/wish/to/place/your/project
git clone git@git.loschmidt.cz/fireprot-asr-vae
cd fireprot-asr-vae
```
Create conda environment for the project
```
conda create --name fireprot-asr-vae python=3.10
```
Activate environment and install requirement and make our project visible in the environment
```
conda activate fireprot-asr-vae
conda install --file requirements.txt
pip install -e .  # by this you make our latents library visible in the environment
```
The framework using PyTorch for Deep Learning therefor we need to install it.
```
conda install pytorch torchvision -c pytorch  # if your work station does not have GPU
```
or with GPU acceleration
```
conda install pytorch torchvision cudatoolkit -c pytorch
```
## Configuration
`Yaml` configuration file is required when running experiments. To make setup of experiments more 
cloud friendly the `cli/config.py` client and `configs/config_template.yaml` exist.

Copy template file:
```
cp configs/config_template.yaml my_cofig.yaml
```
Set all fields by current experiment requirements, e.g. set out dir
```
python path_to_cli/configs.py set out_dir 'experiment1' path_to/my_cofig.yaml
```

## Running

To run the pipeline you have to go to the folder where your experiments will be stored e.g. experiments
```
cd folder/with/experiments
```
### MSA preprocessing
Prepare MSA data to train VAE by running:
```
python path/to/the/main/cli/main.py process /path/to/the/conf.yaml
```
All processing information are logged in the experiment directory folder  `experiment_path_dir/logs/msa_log.txt`

### Training
Then you may modify the configuration file for training options (e.g. custom number of epoch or architecture)
and you can run **TRAINING** by hitting this command to the command line:
```
python path/to/the/main/cli/main.py train /path/to/the/conf.yaml
```
The 80/20% train/test data split is applied, with query sequence left for the training. 
The `batch_size` the size of the dataset from the default but can modify by `conf`.
Training can run for `num_of_epochs` or until the improvement is not reached in last 10 epochs.
Also, linear decrease of L2 regularization factor for 1st quarter (in the case of `num_of_epochs`)
or for first 1000 epochs in unset `num_of_epoch` regime.

### Generating ancestors
Generate ancestors from the latent space is possible via running these commands:
```
python path/to/the/main/cli/main.py ancestors /path/to/the/conf.yaml   # straight evolution protocol
```
Results of generated ancestors with the profile are located in `results/` folder of `experiment_path` directory

## Data transformation

You can transform custom data to the binary:

```
python cli/tasks.py binarize msa_path out_binary_path
```