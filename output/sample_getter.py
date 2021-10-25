import math
import numpy as np


# all models
def get_sampling_frequency(algorithm_name, output_directory, sampling_frequency, temperature_directory):
    if sampling_frequency is None:
        return 1.0 / get_physical_time_step(algorithm_name, output_directory, temperature_directory)
    return sampling_frequency


def get_physical_time_step(algorithm_name, output_directory, temperature_directory):
    acceptance_rates = get_acceptance_rates(output_directory, temperature_directory)
    if algorithm_name == "hxy-gaussian-noise-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis":
        return 2.0 * acceptance_rates[1] * acceptance_rates[0] ** 2  # emergent Langevin diffusivity for Gaussian noise
    return acceptance_rates[1] * acceptance_rates[0] ** 2 / 6.0  # emergent Langevin diffusivity for uniform noise


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
def get_non_normalised_cartesian_magnetisation(output_directory, temperature_directory, beta, no_of_sites):
    return get_entire_sample(output_directory, temperature_directory)[:, 1:3]


def get_magnetisation_norm(output_directory, temperature_directory, beta, no_of_sites):
    return np.linalg.norm(get_non_normalised_cartesian_magnetisation(output_directory, temperature_directory, beta,
                                                                     no_of_sites), axis=1) / no_of_sites


def get_magnetisation_phase(output_directory, temperature_directory, beta, no_of_sites):
    return np.array([math.atan(observation[1] / observation[0]) for observation in
                     get_non_normalised_cartesian_magnetisation(output_directory, temperature_directory, beta,
                                                                no_of_sites)])


def get_cartesian_magnetisation(output_directory, temperature_directory, beta, no_of_sites):
    return get_non_normalised_cartesian_magnetisation(output_directory, temperature_directory, beta,
                                                      no_of_sites) / no_of_sites


def get_magnetic_susceptibility(output_directory, temperature_directory, beta, no_of_sites):
    magnetisation_norm = get_magnetisation_norm(output_directory, temperature_directory, beta, no_of_sites)
    return beta * no_of_sites * (magnetisation_norm - np.mean(magnetisation_norm)) ** 2


def get_non_normalised_total_vortex_polarisation(output_directory, temperature_directory, beta, no_of_sites):
    return get_entire_sample(output_directory, temperature_directory)[:, 3:5]


def get_inverse_vacuum_permittivity(output_directory, temperature_directory, beta, no_of_sites):
    return np.mean(get_entire_sample(output_directory, temperature_directory)[:, 5:7], axis=1) / no_of_sites


def get_helicity_modulus(output_directory, temperature_directory, beta, no_of_sites):
    non_normalised_total_vortex_polarisation = get_non_normalised_total_vortex_polarisation(output_directory,
                                                                                            temperature_directory,
                                                                                            beta, no_of_sites)
    return (get_inverse_vacuum_permittivity(output_directory, temperature_directory, beta, no_of_sites) - beta
            * np.mean((non_normalised_total_vortex_polarisation
                       - np.mean(non_normalised_total_vortex_polarisation, axis=0)) ** 2, axis=1) / no_of_sites)


def get_total_vortex_polarisation(output_directory, temperature_directory, beta, no_of_sites):
    return get_non_normalised_total_vortex_polarisation(output_directory, temperature_directory, beta,
                                                        no_of_sites) / no_of_sites


def get_hxy_topological_sector(output_directory, temperature_directory, beta, no_of_sites):
    non_normalised_total_vortex_polarisation = get_non_normalised_total_vortex_polarisation(output_directory,
                                                                                            temperature_directory,
                                                                                            beta, no_of_sites)
    return compute_hxy_and_maggs_topological_sector(non_normalised_total_vortex_polarisation, no_of_sites)


# Maggs-electrolyte models
def get_sum_of_electric_field(output_directory, temperature_directory, beta, no_of_sites):
    return get_entire_sample(output_directory, temperature_directory)[:, 1:3]


def get_total_polarisation(output_directory, temperature_directory, beta, no_of_sites):
    return get_sum_of_electric_field(output_directory, temperature_directory, beta, no_of_sites) / no_of_sites


def get_inverse_permittivity(output_directory, temperature_directory, beta, no_of_sites):
    sum_of_electric_field_sample = get_sum_of_electric_field(output_directory, temperature_directory, beta, no_of_sites)
    return 1.0 - beta * np.mean((sum_of_electric_field_sample - np.mean(sum_of_electric_field_sample, axis=0)) ** 2,
                                axis=1) / no_of_sites


def get_topological_sector(output_directory, temperature_directory, beta, no_of_sites):
    sum_of_electric_field = get_sum_of_electric_field(output_directory, temperature_directory, beta, no_of_sites)
    return compute_hxy_and_maggs_topological_sector(sum_of_electric_field, no_of_sites)


def get_topological_sector_fluctuations(output_directory, temperature_directory, beta, no_of_sites):
    # n.b. in PRB 91, 155412, topological_sector_fluctuations were defined w/out the factor of 1/2 to account for each
    # Cartesian dimension of the system; in the return line below, this corresponds to the first np.mean() -> np.sum()
    topological_sector_sample = get_topological_sector(output_directory, temperature_directory, beta, no_of_sites)
    return 4.0 * beta * math.pi ** 2 * np.mean((topological_sector_sample - np.mean(topological_sector_sample, axis=0))
                                               ** 2, axis=1)


def compute_hxy_and_maggs_topological_sector(non_normalised_total_polarisation, no_of_sites):
    topological_sector = (non_normalised_total_polarisation + math.pi * no_of_sites ** 0.5) // (2.0 * math.pi
                                                                                                * no_of_sites ** 0.5)
    # correct topological sector so that Ebar_{p, x / y} \in (-pi / L, pi / L] (rather than [-pi / L, pi / L]), where
    # Ebar = Ebar_p + Ebar_w, with Ebar / Ebar_p the total / minimal polarisation on the torus
    length_independent_total_polarisation = non_normalised_total_polarisation / (no_of_sites ** 0.5)
    for index, observation in enumerate(topological_sector):
        for cartesian_index, component in enumerate(observation):
            if component < 0 and (abs(abs(length_independent_total_polarisation[index, cartesian_index])
                                      + component * 3.141592653589793) < 1.0e-12):
                topological_sector[index, cartesian_index] += 1
    return topological_sector
