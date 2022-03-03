import os
import argparse


def is_valid_directory(arg):
    """
    Checks directory path given as input argument exists. Throws an error if it does not.
    """
    if os.path.isdir(arg):
        return arg
    else:
        raise argparse.ArgumentTypeError("{} is not a valid directory.".format(arg))

