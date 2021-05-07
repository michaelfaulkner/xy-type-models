import importlib
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import sys

# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
config_data_getter = importlib.import_module('config_data_getter')
sample_getter = importlib.import_module('sample_getter')
markov_chain_diagnostics = importlib.import_module('markov_chain_diagnostics')
polyspectra = importlib.import_module('polyspectra')


def main(config_file_name, power_spectrum_string):
    basic_config_data = config_data_getter.get_basic_data(config_file_name)
    (algorithm_name, output_directory, integer_lattice_length, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = (basic_config_data[0], basic_config_data[1],
                                                                     basic_config_data[2], basic_config_data[3],
                                                                     basic_config_data[5], basic_config_data[6],
                                                                     basic_config_data[7], basic_config_data[8])

    if ((algorithm_name == 'elementary-electrolyte' or algorithm_name == 'multivalued-electrolyte') and
            (power_spectrum_string == 'magnetisation_norm' or power_spectrum_string == 'helicity_modulus' or
             power_spectrum_string == 'inverse_vacuum_permittivity' or power_spectrum_string ==
             'toroidal_vortex_polarisation')):
        print('ConfigurationError: This is an Maggs-electrolyte model: do not give either magnetisation, '
              'helicity_modulus, inverse_vacuum_permittivity or toroidal_vortex_polarisation as the second positional '
              'argument.')
        exit()
    if ((algorithm_name == 'xy-ecmc' or algorithm_name == 'hxy-ecmc' or algorithm_name == 'xy-metropolis' or
         algorithm_name == 'hxy-metropolis') and (power_spectrum_string == 'inverse_permittivity' or
                                                  power_spectrum_string == 'topological_sector_fluctuations' or
                                                  power_spectrum_string == 'toroidal_polarisation')):
        print('ConfigurationError: This is an XY or HXY model: do not give either inverse_permittivity, '
              'topological_sector_fluctuations or toroidal_polarisation as the second positional argument.')
        exit()

    if no_of_temperature_increments == 0:
        magnitude_of_temperature_increments = 0.0
    else:
        magnitude_of_temperature_increments = (final_temperature -
                                               initial_temperature) / no_of_temperature_increments
    no_of_sites = integer_lattice_length ** 2

    matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
    plt.xlabel(r"frequency, $\omega$", fontsize=10, labelpad=10)
    plt.ylabel(r"$ S_X \left( \omega \right)$", fontsize=10, labelpad=10)
    plt.tick_params(axis='both', which='major', labelsize=10, pad=10)

    temperature = final_temperature
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, no_of_temperature_increments + 1)))
    for i in range(no_of_temperature_increments + 1):
        beta = 1.0 / temperature
        temperature_directory = '/temp_eq_' + f'{temperature:.2f}'
        if no_of_jobs == 1:
            power_spectrum = polyspectra.get_power_spectrum(power_spectrum_string, output_directory,
                                                            temperature_directory, beta, no_of_sites,
                                                            no_of_equilibration_sweeps)
        else:
            number_of_cpus = mp.cpu_count()
            pool = mp.Pool(number_of_cpus)
            power_spectra = pool.starmap(
                polyspectra.get_power_spectrum, [(power_spectrum_string, output_directory + "/run_" +
                                                  str(job_number + 1), temperature_directory, beta, no_of_sites,
                                                  no_of_equilibration_sweeps) for job_number in range(no_of_jobs)])
            pool.close()
            power_spectrum = np.mean(np.array(power_spectra), axis=0)
        plt.plot(power_spectrum[0][0:250], power_spectrum[1][0:250], color=next(colors),
                 label=f'temperature = {temperature:.2f}')
        plt.tight_layout()
        temperature -= magnitude_of_temperature_increments
    legend = plt.legend(loc='upper right', fontsize=10)
    legend.get_frame().set_edgecolor('k')
    legend.get_frame().set_lw(1.5)
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
