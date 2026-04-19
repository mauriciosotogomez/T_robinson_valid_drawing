# import utils functions
import utils
import numpy as np
import time
import pandas as pd
import re
import subprocess

# Arguments
from argparse import ArgumentParser

#Solve LP
from pulp import * # solve LP

#######################
# Read Input  #
#######################

parser = ArgumentParser()
parser.add_argument("-e", "--example_number",
                   dest="n_example", default=0,
                   help="use this example from the list ")

parser.add_argument("-f", "--input_file",
                   dest="input_file", default=None,
                   help="use this .csv file as input")

parser.add_argument("-g", "--generate_random",
                    dest="generate_random_size", default=None,
                    help="generate a random matrix with this size")

parser.add_argument("-m", "--generate_center_matrix",
                    dest="generate_random_center", default=None,
                    help="generate a random center matrix with this size")

parser.add_argument("-v",
                    action="store_true",
                    help="verbose")

parser.add_argument("-d",
                    action="store_true",
                    help="use distance matrix")

parser.add_argument("-c",
                    action="store_true",
                    help="use center matrix")

parser.add_argument("-p",
                    action="store_true",
                    help="use a path instead of centipede")

args = parser.parse_args()
