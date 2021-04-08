import importlib
import os
import sys


# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
output_directory = os.path.abspath(this_directory + "/../")
sys.path.insert(0, output_directory)
config_file = importlib.import_module("config_file")
get_sample = importlib.import_module("get_sample")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")


def main(config_file_name, summary_statistic_string):
    if (summary_statistic_string != 'helicity_modulus' and summary_statistic_string != 'magnetisation' and
            summary_statistic_string != 'magnetisation'):
        IOError('Give one of helicity_modulus, magnetisation or magnetisation as the second positional argument.')
    basic_configuration_data = config_file.get_basic_configuration_data(config_file_name)
    (algorithm_name, sample_directory, lattice_length, no_of_equilibrium_iterations, initial_temperature,
     final_temperature, no_of_temperature_increments) = (basic_configuration_data[0], basic_configuration_data[1],
                                                         basic_configuration_data[2], basic_configuration_data[3],
                                                         basic_configuration_data[5], basic_configuration_data[6],
                                                         basic_configuration_data[7])
    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature - initial_temperature) / no_of_temperature_increments
    no_of_sites = lattice_length

    temperature = initial_temperature
    output_file = open(sample_directory + '/' + summary_statistic_string + '_vs_temperature.dat', 'w')
    output_file.write("temperature".ljust(20) + '\t' + summary_statistic_string.ljust(20) + '\t' +
                      summary_statistic_string + " error".ljust(20) + '\n')
    for i in range(no_of_temperature_increments + 1):
        beta = 1.0 / temperature
        temperature_directory = '/temp_eq_' + str(format(temperature, '.2f'))
        sample = getattr(get_sample, summary_statistic_string)(sample_directory, temperature_directory, beta,
                                                               no_of_sites)[no_of_equilibrium_iterations:]
        sample_mean, sample_error = markov_chain_diagnostics.get_sample_mean_and_error(sample)
        output_file.write("{0:.2e}".format(temperature).ljust(20) + '\t' + "{0:.14e}".format(sample_mean).ljust(20) +
                          '\t' + "{0:.14e}".format(sample_error).ljust(20) + '\n')
        temperature += magnitude_of_temperature_increments
    output_file.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
