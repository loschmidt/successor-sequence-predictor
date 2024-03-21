"""
Client for a successor prediction project
"""
import click

from successors.final_consensus import predict_level3
from successors.majority_voter import majority_voter
from successors.predictor import predict
from successors.parser_handler import RunSetup
from successors.aa_indexes import prepare_indices, prepare_wt_sequence
from successors.index_consensus import predict_level2


class ImproperlyConfigured(Exception):
    """The given configuration is incomplete or otherwise not usable."""
    pass


@click.group()
def cli():
    pass


@cli.command(name='predict')
@click.argument('file_path')
def predict_cli(file_path):
    """ Setup indexes you want to use during analysis """
    print("Predicting successor sequences")
    try:
        run = RunSetup(file_path)
        print(" Target directory: {}".format(run.out_dir))
    except Exception as e:
        print(f" Cannot fetch configuration from {file_path}")
        print(e)
        exit(2)
    predict(run)
    print("  Done")


@cli.command(name='level2')
@click.argument('file_path')
def level2(file_path):
    """ Predict successor sequences based on the given AA indices over all trees"""
    print("Predicting successor from trees per AA index")
    try:
        run = RunSetup(file_path)
        print(" Target directory: {}".format(run.out_dir))
    except Exception as e:
        print(f" Cannot fetch configuration from {file_path}")
        print(e)
        exit(2)
    prepare_wt_sequence(run)
    prepare_indices(run)
    predict_level2(run)
    print("\n  See metric results and generated sequences per AA index in (out_dir_path)/results/level2")
    print("  Done")


@cli.command(name='level3')
@click.argument('file_path')
def level3(file_path):
    """ Predict successor sequences based on the all AA indices """
    print("Predicting successor from all AA index")
    try:
        run = RunSetup(file_path)
        print(" Target directory: {}".format(run.out_dir))
    except Exception as e:
        print(f" Cannot fetch configuration from {file_path}")
        print(e)
        exit(2)
    predict_level3(run)
    print(f"\n  See metric results and generated sequences per AA index in {run.out_dir}/results/level3")
    majority_voter(run)
    print("  Done")


if __name__ == "__main__":
    cli()
