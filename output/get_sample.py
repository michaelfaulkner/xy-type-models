import numpy as np


def helicity_modulus(sample_directory, temperature_directory, beta, no_of_sites):
    sample_mean_1st_derivative_of_potential = np.mean(np.loadtxt(sample_directory + temperature_directory +
                                                                 '/mean_1st_derivative_of_potential_sample.dat',
                                                                 dtype=float, delimiter=','), axis=1)
    sample_mean_2nd_derivative_of_potential = np.mean(np.loadtxt(sample_directory + temperature_directory +
                                                                 '/mean_2nd_derivative_of_potential_sample.dat',
                                                                 dtype=float, delimiter=','), axis=1)
    return sample_mean_2nd_derivative_of_potential - beta * no_of_sites * sample_mean_1st_derivative_of_potential ** 2


def magnetisation(sample_directory, temperature_directory):
    return np.linalg.norm(np.loadtxt(sample_directory + temperature_directory + '/magnetisation_sample.dat',
                                     dtype=float, delimiter=','), axis=1)


def specific_heat(sample_directory, temperature_directory, beta, no_of_sites):
    potential_sample = np.loadtxt(sample_directory + temperature_directory + '/potential_sample.dat', dtype=float,
                                  delimiter=',')
    return beta ** 2 * (potential_sample - np.mean(potential_sample)) ** 2 / no_of_sites
