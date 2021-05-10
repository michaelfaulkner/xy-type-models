import math
import numpy as np


# all models
def get_acceptance_rates(output_directory, temperature_directory):
    return np.atleast_1d(np.loadtxt(f"{output_directory}/{temperature_directory}/acceptance_rates.csv", dtype=float,
                                    delimiter=","))


def get_no_of_events(output_directory, temperature_directory):
    return np.atleast_1d(np.loadtxt(f"{output_directory}/{temperature_directory}/no_of_events.csv", dtype=float,
                                    delimiter=","))


def get_entire_sample(output_directory, temperature_directory):
    return np.loadtxt(f"{output_directory}/{temperature_directory}/sample.csv", dtype=float, delimiter=",")


def get_potential(output_directory, temperature_directory, beta, no_of_sites):
    return get_entire_sample(output_directory, temperature_directory)[:, 0]


def get_specific_heat(output_directory, temperature_directory, beta, no_of_sites):
    potential_sample = get_potential(output_directory, temperature_directory, beta, no_of_sites)
    return beta ** 2 * (potential_sample - np.mean(potential_sample)) ** 2 / no_of_sites


# xy models
def get_helicity_modulus(output_directory, temperature_directory, beta, no_of_sites):
    sum_of_1st_derivative_of_potential_sample = get_entire_sample(output_directory, temperature_directory)[:, 3:5]
    return (get_inverse_vacuum_permittivity(output_directory, temperature_directory, beta, no_of_sites) - beta
            * np.mean((sum_of_1st_derivative_of_potential_sample - np.mean(sum_of_1st_derivative_of_potential_sample,
                                                                           axis=0)) ** 2, axis=1) / no_of_sites)


def get_inverse_vacuum_permittivity(output_directory, temperature_directory, beta, no_of_sites):
    return np.mean(get_entire_sample(output_directory, temperature_directory)[:, 5:7], axis=1) / no_of_sites


def get_toroidal_vortex_polarisation(output_directory, temperature_directory, beta, no_of_sites):
    return get_entire_sample(output_directory, temperature_directory)[:, 3:5] / no_of_sites


def get_magnetisation_norm(output_directory, temperature_directory, beta, no_of_sites):
    return np.linalg.norm(get_non_normalised_cartesian_magnetisation(output_directory, temperature_directory, beta,
                                                                     no_of_sites), axis=1) / no_of_sites


def get_magnetisation_phase(output_directory, temperature_directory, beta, no_of_sites):
    return np.array([math.atan(observation[1] / observation[0]) for observation in
                     get_non_normalised_cartesian_magnetisation(output_directory, temperature_directory, beta,
                                                                no_of_sites)]) / no_of_sites


def get_non_normalised_cartesian_magnetisation(output_directory, temperature_directory, beta, no_of_sites):
    return get_entire_sample(output_directory, temperature_directory)[:, 1:3]


# Maggs-electrolyte models
def get_sum_of_electric_field(output_directory, temperature_directory, beta, no_of_sites):
    return get_entire_sample(output_directory, temperature_directory)[:, 1:3]


def get_toroidal_polarisation(output_directory, temperature_directory, beta, no_of_sites):
    return get_sum_of_electric_field(output_directory, temperature_directory, beta,
                                     no_of_sites) / no_of_sites


def get_inverse_permittivity(output_directory, temperature_directory, beta, no_of_sites):
    sum_of_electric_field_sample = get_sum_of_electric_field(output_directory, temperature_directory, beta,
                                                             no_of_sites)
    return 1.0 - beta * np.mean((sum_of_electric_field_sample - np.mean(sum_of_electric_field_sample, axis=0)) ** 2,
                                axis=1) / no_of_sites


def get_topological_sector_fluctuations(output_directory, temperature_directory, beta, no_of_sites):
    # n.b. in PRB 91, 155412, topological_sector_fluctuations were defined w/out the factor of 1/2 to account for each
    # Cartesian dimension of the system; in the return line below, this corresponds to the first np.mean() -> np.sum()
    topological_sector_sample = get_topological_sector(output_directory, temperature_directory, beta,
                                                       no_of_sites)
    return 4.0 * beta * math.pi ** 2 * np.mean((topological_sector_sample - np.mean(topological_sector_sample, axis=0))
                                               ** 2, axis=1)


def get_topological_sector(output_directory, temperature_directory, beta, no_of_sites):
    sum_of_electric_field_sample = get_sum_of_electric_field(output_directory, temperature_directory, beta,
                                                             no_of_sites)
    return (sum_of_electric_field_sample + math.pi * no_of_sites ** 0.5) // (2.0 * math.pi * no_of_sites ** 0.5)
