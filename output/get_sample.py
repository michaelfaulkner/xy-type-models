import numpy as np


def acceptance_rates(sample_directory, temperature_directory):
    return np.atleast_1d(np.loadtxt(sample_directory + temperature_directory + '/acceptance_rates.dat', dtype=float,
                                    delimiter=','))


def number_of_events(sample_directory, temperature_directory):
    return np.atleast_1d(np.loadtxt(sample_directory + temperature_directory + '/number_of_events.dat', dtype=float,
                                    delimiter=','))


def inverse_permittivity(sample_directory, temperature_directory, beta, number_of_sites):
    harmonic_mode_sample = np.mean(np.loadtxt(sample_directory + temperature_directory +
                                              '/sum_of_electric_field_x_sample.dat', dtype=float, delimiter=','),
                                   axis=1)
    return 1.0 - beta * (harmonic_mode_sample - np.mean(harmonic_mode_sample)) ** 2 / number_of_sites


def helicity_modulus(sample_directory, temperature_directory, beta, number_of_sites):
    mean_1st_derivative_of_potential_sample = np.mean(np.loadtxt(sample_directory + temperature_directory +
                                                                 '/mean_1st_derivative_of_potential_sample.dat',
                                                                 dtype=float, delimiter=','), axis=1)
    mean_2nd_derivative_of_potential_sample = np.mean(np.loadtxt(sample_directory + temperature_directory +
                                                                 '/mean_2nd_derivative_of_potential_sample.dat',
                                                                 dtype=float, delimiter=','), axis=1)
    return mean_2nd_derivative_of_potential_sample - beta * (mean_1st_derivative_of_potential_sample - np.mean(
        mean_1st_derivative_of_potential_sample)) ** 2 / number_of_sites


def magnetisation(sample_directory, temperature_directory, beta, number_of_sites):
    return np.linalg.norm(np.loadtxt(sample_directory + temperature_directory + '/magnetisation_sample.dat',
                                     dtype=float, delimiter=','), axis=1) / number_of_sites


def specific_heat(sample_directory, temperature_directory, beta, number_of_sites):
    potential_sample = np.loadtxt(sample_directory + temperature_directory + '/potential_sample.dat', dtype=float,
                                  delimiter=',')
    return beta ** 2 * (potential_sample - np.mean(potential_sample)) ** 2 / number_of_sites
