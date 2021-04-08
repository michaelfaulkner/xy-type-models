import csv
import importlib
import numpy as np
import os
import sys

config_file = importlib.import_module("config_file")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")


def main(config_file_name):
    (algorithm_name, sample_directory, lattice_length, no_of_equilibrium_iterations, no_of_observations,
     initial_temperature, final_temperature, no_of_temperature_increments) = config_file.get_basic_configuration_data(config_file_name)


if __name__ == '__main__':
    main(sys.argv[1])
