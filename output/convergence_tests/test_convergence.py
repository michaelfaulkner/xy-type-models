import csv
import importlib
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys


def main(algorithm_name, config_file_name):
    matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
    with open(config_file_name, 'r') as config_file:
        config_data = csv.reader(config_file, delimiter='\t')
        for index, row in enumerate(config_data):
            if index == 0:
                sample_directory = row[0].replace("'", "").replace("output directory", "").replace(" ", "")
            if index == 2:
                no_of_equilibrium_iterations = int(row[0].replace("No. of equilibration iterations", "").replace(" ", ""))
            if index == 4:
                temperature = float(row[0].replace("Initial temperature", "").replace(" ", ""))
                temperature_directory = '/temp_eq_' + str(format(temperature, '.2f'))
                break

    if algorithm_name == 'elementary-electrolyte' or algorithm_name == 'multivalued-electrolyte':
        sample = np.loadtxt(sample_directory + temperature_directory + '/potential_sample.dat', dtype=float,
                            delimiter=',')
    elif (algorithm_name == 'hxy-ecmc' or algorithm_name == 'hxy-metropolis' or algorithm_name == 'xy-ecmc' or
            algorithm_name == 'xy-metropolis'):
        magnetisation_x = np.loadtxt(sample_directory + temperature_directory + '/magnetisation_x_sample.dat',
                                     dtype=float)[no_of_equilibrium_iterations:]
        magnetisation_y = np.loadtxt(sample_directory + temperature_directory + '/magnetisation_y_sample.dat',
                                     dtype=float)[no_of_equilibrium_iterations:]
        sample = (magnetisation_x ** 2 + magnetisation_y ** 2) ** 0.5
    else:
        IOError('Give one of elementary-electrolyte, multivalued-electrolyte, hxy-ecmc, hxy-metropolis, xy-ecmc or '
                'xy-metropolis as the first positional argument of the test_convergence.py script.')

    if algorithm_name == 'elementary-electrolyte':
        reference_sample = np.loadtxt('output/convergence_tests/maggs_electrolyte/elementary_charges/elementary_maggs_'
                                      'electrolyte_8x8_square_temp_1_point_5_potential_reference_sample.dat',
                                      dtype=float)
    if algorithm_name == 'multivalued-electrolyte':
        reference_sample = np.loadtxt('output/convergence_tests/maggs_electrolyte/multivalued_charges/multivalued_'
                                      'maggs_electrolyte_8x8_square_temp_1_point_5_potential_reference_sample.dat',
                                      dtype=float)
    if algorithm_name == 'hxy-ecmc' or algorithm_name == 'hxy-metropolis':
        reference_sample = np.loadtxt('output/convergence_tests/hxy/hxy_8x8_sites_temp_1_point_3_magnetisation_norm_'
                                      'reference_sample.dat', dtype=float)
    if algorithm_name == 'xy-ecmc' or algorithm_name == 'xy-metropolis':
        reference_sample = np.loadtxt('output/convergence_tests/xy/xy_8x8_sites_temp_0_point_8_magnetisation_norm_'
                                      'reference_sample.dat', dtype=float)

    # Add the directory that contains the module plotting_functions to sys.path
    this_directory = os.path.dirname(os.path.abspath(__file__))
    output_directory = os.path.abspath(this_directory + "/../")
    sys.path.insert(0, output_directory)
    markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")
    effective_sample_size = markov_chain_diagnostics.get_effective_sample_size(sample)
    print(f"Effective sample size = {effective_sample_size} (from a total sample size of {len(sample)}).")

    '''magnetisation_2_x = np.loadtxt('output/convergence_tests/hxy/metropolis/temp_eq_1.50/magn_x_sample.dat',
                                   dtype=float)
    magnetisation_2_y = np.loadtxt('output/convergence_tests/hxy/metropolis/temp_eq_1.50/magn_y_sample.dat',
                                   dtype=float)
    sample_2 = (magnetisation_2_x ** 2 + magnetisation_2_y ** 2) ** 0.5
    sample_2_cdf = get_cumulative_distribution(sample_2)
    plt.plot(sample_2_cdf[0], sample_2_cdf[1], color='r', linewidth=6, linestyle='-', label='metropolis data')'''

    reference_cdf = get_cumulative_distribution(reference_sample)
    sample_cdf = get_cumulative_distribution(sample)

    plt.plot(reference_cdf[0], reference_cdf[1], color='r', linewidth=4, linestyle='-', label='reference data')
    plt.plot(sample_cdf[0], sample_cdf[1], color='k', linewidth=2, linestyle='-', label='xy-type-models data')

    plt.xlabel(r"$x$", fontsize=15, labelpad=10)
    plt.ylabel(r"$ \mathbb{P} \left( X < x \right)$", fontsize=15, labelpad=10)
    plt.tick_params(axis='both', which='major', labelsize=14, pad=10)
    legend = plt.legend(loc='upper left', fontsize=10)
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_lw(1.5)
    plt.tight_layout()
    plt.show()


def get_cumulative_distribution(input_sample):
    bin_values = np.arange(1, len(input_sample) + 1) / float(len(input_sample))
    ordered_input_sample = np.sort(input_sample)
    return ordered_input_sample, bin_values


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
