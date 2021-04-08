import csv
import importlib
import numpy as np
import os
import sys


# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
output_directory = os.path.abspath(this_directory + "/../")
sys.path.insert(0, output_directory)
config_file = importlib.import_module("config_file")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")


def main(config_file_name):
    basic_configuration_data = config_file.get_basic_configuration_data(config_file_name)
    (algorithm_name, sample_directory, no_of_equilibrium_iterations, initial_temperature, final_temperature,
     no_of_temperature_increments) = (basic_configuration_data[0], basic_configuration_data[1],
                                      basic_configuration_data[3], basic_configuration_data[5],
                                      basic_configuration_data[6], basic_configuration_data[7])


if __name__ == '__main__':
    main(sys.argv[1])
