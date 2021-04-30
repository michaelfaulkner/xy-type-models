import math
import numpy as np


# all models
def get_acceptance_rates(simulation_directory, temperature_directory):
    return np.atleast_1d(np.loadtxt(simulation_directory + temperature_directory + '/acceptance_rates.csv', dtype=float,
                                    delimiter=','))


def get_number_of_events(simulation_directory, temperature_directory):
    return np.atleast_1d(np.loadtxt(simulation_directory + temperature_directory + '/number_of_events.csv', dtype=float,
                                    delimiter=','))


def get_entire_sample(simulation_directory, temperature_directory):
    return np.loadtxt(simulation_directory + temperature_directory + '/sample.csv', dtype=float, delimiter=',')


def get_potential(simulation_directory, temperature_directory, beta, number_of_sites):
    return get_entire_sample(simulation_directory, temperature_directory)[:, 0]


def get_specific_heat(simulation_directory, temperature_directory, beta, number_of_sites):
    potential_sample = get_potential(simulation_directory, temperature_directory, beta, number_of_sites)
    return beta ** 2 * (potential_sample - np.mean(potential_sample)) ** 2 / number_of_sites


# xy models
def get_helicity_modulus(simulation_directory, temperature_directory, beta, number_of_sites):
    sum_of_1st_derivative_of_potential_sample = get_entire_sample(simulation_directory, temperature_directory)[:, 3:5]
    return (get_inverse_vacuum_permittivity(simulation_directory, temperature_directory, beta, number_of_sites) - beta
            * np.mean((sum_of_1st_derivative_of_potential_sample - np.mean(sum_of_1st_derivative_of_potential_sample,
                                                                           axis=0)) ** 2, axis=1) / number_of_sites)


def get_inverse_vacuum_permittivity(simulation_directory, temperature_directory, beta, number_of_sites):
    return np.mean(get_entire_sample(simulation_directory, temperature_directory)[:, 5:7], axis=1) / number_of_sites


def get_toroidal_vortex_polarisation(simulation_directory, temperature_directory, beta, number_of_sites):
    return get_entire_sample(simulation_directory, temperature_directory)[:, 3:5] / number_of_sites


def get_magnetisation_norm(simulation_directory, temperature_directory, beta, number_of_sites):
    return np.linalg.norm(get_non_normalised_cartesian_magnetisation(simulation_directory, temperature_directory, beta,
                                                                     number_of_sites), axis=1) / number_of_sites


def get_non_normalised_cartesian_magnetisation(simulation_directory, temperature_directory, beta, number_of_sites):
    return get_entire_sample(simulation_directory, temperature_directory)[:, 1:3]


# Maggs-electrolyte models
def get_sum_of_electric_field(simulation_directory, temperature_directory, beta, number_of_sites):
    return get_entire_sample(simulation_directory, temperature_directory)[:, 1:3]


def get_toroidal_polarisation(simulation_directory, temperature_directory, beta, number_of_sites):
    return get_sum_of_electric_field(simulation_directory, temperature_directory, beta,
                                     number_of_sites) / number_of_sites


def get_inverse_permittivity(simulation_directory, temperature_directory, beta, number_of_sites):
    sum_of_electric_field_sample = get_sum_of_electric_field(simulation_directory, temperature_directory, beta,
                                                             number_of_sites)
    return 1.0 - beta * np.mean((sum_of_electric_field_sample - np.mean(sum_of_electric_field_sample, axis=0)) ** 2,
                                axis=1) / number_of_sites


def get_topological_sector_fluctuations(simulation_directory, temperature_directory, beta, number_of_sites):
    # n.b. in PRB 91, 155412, topological_sector_fluctuations were defined w/out the factor of 1/2 to account for each
    # Cartesian dimension of the system; in the return line below, this corresponds to the first np.mean() -> np.sum()
    topological_sector_sample = get_topological_sector(simulation_directory, temperature_directory, beta,
                                                       number_of_sites)
    return 4.0 * beta * math.pi ** 2 * np.mean((topological_sector_sample - np.mean(topological_sector_sample, axis=0))
                                               ** 2, axis=1)


def get_topological_sector(simulation_directory, temperature_directory, beta, number_of_sites):
    sum_of_electric_field_sample = get_sum_of_electric_field(simulation_directory, temperature_directory, beta,
                                                             number_of_sites)
    return (sum_of_electric_field_sample + math.pi * number_of_sites ** 0.5) // (2.0 * math.pi * number_of_sites ** 0.5)
