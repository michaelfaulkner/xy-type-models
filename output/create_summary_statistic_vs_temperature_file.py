import importlib
import os
import sys


# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
output_directory = os.path.abspath(this_directory + '/../')
sys.path.insert(0, output_directory)
config_file = importlib.import_module('config_file')
get_sample = importlib.import_module('get_sample')
markov_chain_diagnostics = importlib.import_module('markov_chain_diagnostics')


def main(config_file_name, summary_statistic_string):
    if (summary_statistic_string != 'acceptance_rates' and summary_statistic_string != 'number_of_events' and
            summary_statistic_string != 'helicity_modulus' and summary_statistic_string != 'magnetisation' and
            summary_statistic_string != 'specific_heat'):
        print('Give one of acceptance_rates, number_of_events, helicity_modulus, magnetisation or specific_heat as the '
              'second positional argument.')
        exit()

    basic_configuration_data = config_file.get_basic_configuration_data(config_file_name)
    (algorithm_name, sample_directory, lattice_length, number_of_equilibrium_iterations, initial_temperature,
     final_temperature, number_of_temperature_increments) = (basic_configuration_data[0], basic_configuration_data[1],
                                                             basic_configuration_data[2], basic_configuration_data[3],
                                                             basic_configuration_data[5], basic_configuration_data[6],
                                                             basic_configuration_data[7])

    if (algorithm_name == 'xy-ecmc' or algorithm_name == 'hxy-ecmc') and summary_statistic_string == 'acceptance_rates':
        print('ConfigurationError: These samples were generated by an event-chain Monte Carlo algorithm: rather than '
              'acceptance_rates, give number_of_events as the second positional argument.')
        exit()
    if ((algorithm_name == 'xy-metropolis' or algorithm_name == 'hxy-metropolis' or
         algorithm_name == 'elementary-electrolyte' or algorithm_name == 'multivalued-electrolyte') and
            summary_statistic_string == 'number_of_events'):
        print('ConfigurationError: These samples were generated by a Metropolis algorithm: rather than '
              'number_of_events, give acceptance_rates as the second positional argument.')
        exit()
    if ((algorithm_name == 'elementary-electrolyte' or algorithm_name == 'multivalued-electrolyte') and
            (summary_statistic_string == 'magnetisation' or summary_statistic_string == 'helicity_modulus')):
        print('ConfigurationError: This is an Maggs-electrolyte model: do not give either magnetisation or '
              'helicity_modulus as the second positional argument.')
        exit()

    if number_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature -
                                               initial_temperature) / number_of_temperature_increments
    number_of_sites = lattice_length ** 2

    temperature = initial_temperature
    output_file = open(sample_directory + '/' + summary_statistic_string + '_vs_temperature.dat', 'w')
    if summary_statistic_string == 'acceptance_rates':
        temperature_directory = '/temp_eq_' + '{0:.2f}'.format(temperature)
        acceptance_rates = get_sample.acceptance_rates(sample_directory, temperature_directory)
        if len(acceptance_rates) == 1:
            output_file.write('temperature'.ljust(15) + 'rotational acceptance rate' + '\n')
        elif len(acceptance_rates) == 2:
            if algorithm_name == 'elementary-electrolyte' or algorithm_name == 'multivalued-electrolyte':
                output_file.write('temperature'.ljust(15) + 'acceptance rate (rotational moves)'.ljust(40) +
                                  'acceptance rate (charge hops)' + '\n')
            else:
                output_file.write('temperature'.ljust(15) + 'acceptance rate (rotational moves)'.ljust(40) +
                                  'acceptance rate (external global moves)' + '\n')
        else:
            output_file.write('temperature'.ljust(15) + 'acceptance rate (field rotations)'.ljust(40) +
                              'acceptance rate (charge hops)'.ljust(40) + 'acceptance rate (global moves)' + '\n')
    elif summary_statistic_string == 'number_of_events':
        temperature_directory = '/temp_eq_' + '{0:.2f}'.format(temperature)
        number_of_events = get_sample.number_of_events(sample_directory, temperature_directory)
        if len(number_of_events) == 1:
            output_file.write('temperature'.ljust(15) + 'number of events (field rotations)' + '\n')
        else:
            output_file.write('temperature'.ljust(15) + 'number of events (field rotations)'.ljust(40) +
                              'acceptance rate (external global moves)' + '\n')
    else:
        output_file.write('temperature'.ljust(15) + summary_statistic_string.ljust(25) + summary_statistic_string +
                          ' error' + '\n')

    for i in range(number_of_temperature_increments + 1):
        beta = 1.0 / temperature
        temperature_directory = '/temp_eq_' + '{0:.2f}'.format(temperature)
        if summary_statistic_string == 'acceptance_rates' or summary_statistic_string == 'number_of_events':
            get_sample_method = getattr(get_sample, summary_statistic_string)
            acceptance_rates_or_number_of_events = get_sample_method(sample_directory, temperature_directory)
            if len(acceptance_rates_or_number_of_events) == 1:
                output_file.write('{0:.2e}'.format(temperature).ljust(15) +
                                  '{0:.14e}'.format(acceptance_rates_or_number_of_events[0]) + '\n')
            elif len(acceptance_rates_or_number_of_events) == 2:
                output_file.write('{0:.2e}'.format(temperature).ljust(15) +
                                  '{0:.14e}'.format(acceptance_rates_or_number_of_events[0]).ljust(40) +
                                  '{0:.14e}'.format(acceptance_rates_or_number_of_events[1]) + '\n')
            else:
                output_file.write('{0:.2e}'.format(temperature).ljust(15) +
                                  '{0:.14e}'.format(acceptance_rates_or_number_of_events[0]).ljust(40) +
                                  '{0:.14e}'.format(acceptance_rates_or_number_of_events[1]).ljust(40) +
                                  '{0:.14e}'.format(acceptance_rates_or_number_of_events[2]) + '\n')
        else:
            get_sample_method = getattr(get_sample, summary_statistic_string)
            sample = get_sample_method(sample_directory, temperature_directory, beta, number_of_sites)[
                     number_of_equilibrium_iterations:]
            sample_mean, sample_error = markov_chain_diagnostics.get_sample_mean_and_error(sample)
            output_file.write('{0:.2e}'.format(temperature).ljust(15) + '{0:.14e}'.format(sample_mean).ljust(25) +
                              '{0:.14e}'.format(sample_error) + '\n')
        temperature += magnitude_of_temperature_increments
    output_file.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
