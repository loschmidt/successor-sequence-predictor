__author__ = "Pavel Kohout <xkohou15@vutbr.cz>"
__date__ = "2024/02/16"
__description__ = "Importing the and parsing of configuration file and setup run environment"

import datetime
import os

import yaml


class RunSetup:
    def __init__(self, config_file=None):
        """
        Parse and prepare structure for prediction of successor
        :param config_file: path to the yaml configuration path
        """
        with open(config_file, "r") as f:
            conf = yaml.safe_load(f)

        # sanity check of the keys in the configuration file
        for k in {'out_dir', 'input', 'validation', 'transition'}:  # out_dir is always present
            if k not in conf.keys():
                print(f"  The required key {k} is not used in configuration file {config_file}")
                exit(1)

        # setup directory structure if needed
        directories = ["tmp_files", "results", "config", "logs", "index_fld", "ground_truth", "fasta"]
        for exp_dir in directories:
            attr_dir = os.path.join(conf['out_dir'], exp_dir)
            os.makedirs(attr_dir, exist_ok=True)
            setattr(self, exp_dir, attr_dir)
        # set the rest of the keys
        for k, v in conf.items():
            setattr(self, k, v)
        # copy config file to the results directory
        now = datetime.datetime.now()
        file_name = now.strftime("%Y-%m-%d_%H-%M.json")
        yaml.dump(conf, open(os.path.join(self.config, f"{file_name}"), 'w'))
        print(f"Configuration file stored in {os.path.join(self.config, file_name)}")

    def __getattr__(self, attr_name):
        """Set the default values if not specified in the configuration file"""
        fallback_values = {
            "protein_name": "protein_name",
            "query": "query",
            "highlight_top": 20,
            "highlight_pos": [],
            "conservation_offset": None
        }
        if attr_name in fallback_values:
            return fallback_values[attr_name]
        else:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{attr_name}'")
