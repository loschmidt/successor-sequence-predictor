"""
Client for a successor prediction project
"""
import click
from successors.predictor import predict
from successors.parser_handler import RunSetup


class ImproperlyConfigured(Exception):
    """The given configuration is incomplete or otherwise not usable."""
    pass


@click.group()
def cli():
    pass


@cli.command(name='predict')
@click.argument('file_path')
def configure_indexes(file_path):
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


if __name__ == "__main__":
    cli()
