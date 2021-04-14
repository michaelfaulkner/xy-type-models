import math
import numpy as np


# all models
def acceptance_rates(sample_directory, temperature_directory):
    return np.atleast_1d(np.loadtxt(sample_directory + temperature_directory + '/acceptance_rates.dat', dtype=float,
                                    delimiter=','))


def number_of_events(sample_directory, temperature_directory):
    return np.atleast_1d(np.loadtxt(sample_directory + temperature_directory + '/number_of_events.dat', dtype=float,
                                    delimiter=','))


def specific_heat(sample_directory, temperature_directory, beta, number_of_sites):
    potential_sample = np.loadtxt(sample_directory + temperature_directory + '/potential_sample.dat', dtype=float,
                                  delimiter=',')
    return beta ** 2 * (potential_sample - np.mean(potential_sample)) ** 2 / number_of_sites


# xy models
def helicity_modulus(sample_directory, temperature_directory, beta, number_of_sites):
    sum_of_1st_derivative_of_potential_sample = np.loadtxt(sample_directory + temperature_directory +
                                                           '/sum_of_1st_deriv_of_potential_sample.dat', dtype=float,
                                                           delimiter=',')
    sum_of_2nd_derivative_of_potential_sample = np.loadtxt(sample_directory + temperature_directory +
                                                           '/sum_of_2nd_deriv_of_potential_sample.dat', dtype=float,
                                                           delimiter=',')
    return np.mean(sum_of_2nd_derivative_of_potential_sample, axis=1) - beta * np.mean(
        (sum_of_1st_derivative_of_potential_sample - np.mean(sum_of_1st_derivative_of_potential_sample, axis=0)) ** 2,
        axis=1) / number_of_sites


def magnetisation(sample_directory, temperature_directory, beta, number_of_sites):
    return np.linalg.norm(np.loadtxt(sample_directory + temperature_directory + '/magnetisation_sample.dat',
                                     dtype=float, delimiter=','), axis=1) / number_of_sites


# Maggs-electrolyte models
def inverse_permittivity(sample_directory, temperature_directory, beta, number_of_sites):
    sum_of_electric_field_sample = np.loadtxt(sample_directory + temperature_directory + '/field_sum_sample.dat',
                                              dtype=float, delimiter=',')
    return 1.0 - beta * np.mean((sum_of_electric_field_sample - np.mean(sum_of_electric_field_sample, axis=0)) ** 2,
                                axis=1) / number_of_sites


def topological_sector_fluctuations(sample_directory, temperature_directory, beta, number_of_sites):
    topological_sector_sample = topological_sector(sample_directory, temperature_directory, beta, number_of_sites)
    return 4.0 * beta * math.pi ** 2 * np.mean(
        (topological_sector_sample - np.mean(topological_sector_sample, axis=0)) ** 2, axis=1) / number_of_sites


def topological_sector(sample_directory, temperature_directory, beta, number_of_sites):
    sum_of_electric_field_sample = np.loadtxt(sample_directory + temperature_directory + '/field_sum_sample.dat',
                                              dtype=float, delimiter=',')
    return (sum_of_electric_field_sample + math.pi * number_of_sites ** 0.5) // (2.0 * math.pi * number_of_sites ** 0.5)
