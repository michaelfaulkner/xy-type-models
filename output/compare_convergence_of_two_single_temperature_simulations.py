import importlib
import matplotlib
import matplotlib.pyplot as plt
import os
import sys


# Add the directory that contains config_file, sample_getter and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
config_file = importlib.import_module("config_file")
sample_getter = importlib.import_module("sample_getter")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")


def main(config_file_name_1, config_file_name_2):
    matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

    basic_configuration_data_1 = config_file.get_basic_configuration_data(config_file_name_1)
    (algorithm_name_1, simulation_directory_1, lattice_length_1, no_of_equilibrium_iterations_1, temperature_1,
     no_of_temperature_increments_1) = (
        basic_configuration_data_1[0], basic_configuration_data_1[1], basic_configuration_data_1[2],
        basic_configuration_data_1[3], basic_configuration_data_1[5], basic_configuration_data_1[7])
    if no_of_temperature_increments_1 != 0:
        print('ConfigurationError: in the first configuration file, the value of the number of temperature increments '
              'does not equal 0. In order to compare to single-temperature simulations, this is required.')
        exit()

    basic_configuration_data_2 = config_file.get_basic_configuration_data(config_file_name_2)
    (algorithm_name_2, simulation_directory_2, lattice_length_2, no_of_equilibrium_iterations_2, temperature_2,
     no_of_temperature_increments_2) = (
        basic_configuration_data_2[0], basic_configuration_data_2[1], basic_configuration_data_2[2],
        basic_configuration_data_2[3], basic_configuration_data_2[5], basic_configuration_data_2[7])
    if no_of_temperature_increments_2 != 0:
        print('ConfigurationError: in the second configuration file, the value of the number of temperature increments '
              'does not equal 0. In order to compare to single-temperature simulations, this is required.')
        exit()

    if (((algorithm_name_1 == 'hxy-ecmc' or algorithm_name_1 == 'hxy-metropolis') and not
            (algorithm_name_2 == 'hxy-ecmc' or algorithm_name_2 == 'hxy-metropolis')) or
        ((algorithm_name_1 == 'xy-ecmc' or algorithm_name_1 == 'xy-metropolis') and not
            (algorithm_name_2 == 'xy-ecmc' or algorithm_name_2 == 'xy-metropolis')) or
        (algorithm_name_1 == 'elementary-electrolyte' and algorithm_name_2 != 'elementary-electrolyte') or
            (algorithm_name_1 == 'multivalued-electrolyte' and algorithm_name_2 != 'multivalued-electrolyte')):
        print('ConfigurationError: give the same model in each configuration file.')
        exit()
    if temperature_1 != temperature_2:
        print('ConfigurationError: give the same initial temperature in each configuration file.')
        exit()
    if lattice_length_1 != lattice_length_2:
        print('ConfigurationError: give the same lattice length in each configuration file.')
        exit()

    beta = 1.0 / temperature_1
    number_of_sites = lattice_length_1 ** 2
    temperature_directory = '/temp_eq_' + str(format(temperature_1, '.2f'))

    if algorithm_name_1 == 'elementary-electrolyte' or algorithm_name_1 == 'multivalued-electrolyte':
        sample_1 = sample_getter.get_potential(simulation_directory_1, temperature_directory, beta, number_of_sites)[
                 no_of_equilibrium_iterations_1:]
        sample_2 = sample_getter.get_potential(simulation_directory_2, temperature_directory, beta, number_of_sites)[
                   no_of_equilibrium_iterations_2:]
    elif (algorithm_name_1 == 'hxy-ecmc' or algorithm_name_1 == 'hxy-metropolis' or algorithm_name_1 == 'xy-ecmc' or
            algorithm_name_1 == 'xy-metropolis'):
        sample_1 = sample_getter.get_magnetisation_norm(simulation_directory_1, temperature_directory, beta, number_of_sites)[
                 no_of_equilibrium_iterations_1:]
        sample_2 = sample_getter.get_magnetisation_norm(simulation_directory_2, temperature_directory, beta,
                                                        number_of_sites)[
                   no_of_equilibrium_iterations_2:]

    effective_sample_size_1 = markov_chain_diagnostics.get_effective_sample_size(sample_1)
    print(f"Effective sample size (first config file) = {effective_sample_size_1} (from a total sample size of "
          f"{len(sample_1)}).")
    effective_sample_size_2 = markov_chain_diagnostics.get_effective_sample_size(sample_2)
    print(f"Effective sample size (second config file) = {effective_sample_size_2} (from a total sample size of "
          f"{len(sample_2)}).")

    sample_1_cdf = markov_chain_diagnostics.get_cumulative_distribution(sample_1)
    sample_2_cdf = markov_chain_diagnostics.get_cumulative_distribution(sample_2)
    plt.plot(sample_1_cdf[0], sample_1_cdf[1], color='r', linewidth=4, linestyle='-', label='first config file')
    plt.plot(sample_2_cdf[0], sample_2_cdf[1], color='k', linewidth=2, linestyle='-', label='second config file')

    plt.xlabel(r"$x$", fontsize=15, labelpad=10)
    plt.ylabel(r"$ \mathbb{P} \left( X < x \right)$", fontsize=15, labelpad=10)
    plt.tick_params(axis='both', which='major', labelsize=14, pad=10)
    legend = plt.legend(loc='upper left', fontsize=10)
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_lw(1.5)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
