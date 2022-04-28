"""This module contains methods that compute the Cramér-von Mises statistics of the magnetisation phase."""
from scipy import stats
import importlib
import math
import numpy as np

sample_getter = importlib.import_module("sample_getter")


def get_cvm_mag_phase_vs_temperature(algorithm_name, output_directory, sample_directory, no_of_equilibration_sweeps,
                                     no_of_observations, temperatures, external_global_moves_string, no_of_jobs, pool,
                                     length):
    try:
        with open(f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}"
                  f"_{external_global_moves_string}_{length}x{length}_sites_{no_of_observations}_obs.tsv",
                  "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            cvms, cvm_errors = output_file_sans_header[1], output_file_sans_header[2]
    except IOError:
        cvm_file = open(
            f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}_"
            f"{external_global_moves_string}_{length}x{length}_sites_{no_of_observations}_obs.tsv", "w")
        cvm_file.write("# temperature".ljust(30) + "n omega_n^2".ljust(30) + "n omega_n^2 error" + "\n")
        cvms, cvm_errors = [], []
        for temperature_index, temperature in enumerate(temperatures):
            print(f"Temperature = {temperature:.4f}")
            cvm_vs_job = pool.starmap(get_cramervonmises_of_magnetisation_phase, [
                (f"{sample_directory}/job_{job_index}", temperature, temperature_index, length ** 2,
                 no_of_equilibration_sweeps) for job_index in range(no_of_jobs)])
            cvm, cvm_error = np.mean(cvm_vs_job), np.std(cvm_vs_job) / len(cvm_vs_job) ** 0.5
            cvms.append(cvm)
            cvm_errors.append(cvm_error)
            cvm_file.write(f"{temperature:.14e}".ljust(30) + f"{cvm:.14e}".ljust(30) + f"{cvm_error:.14e}" + "\n")
        cvm_file.close()
    return cvms, cvm_errors


def get_cramervonmises_of_magnetisation_phase(sample_directory, temperature, temperature_index, no_of_sites,
                                              no_of_equilibration_sweeps):
    r"""
    Returns (for the sample of the magnetisation phase) the normalised Cramér-von Mises statistic
    n \omega_n^2 := n \int (F_n(x) - F(x))^2 dF(x), where n is the sample size and F(x) is the CDF of the continuous
    uniform distribution Uniform(-pi, pi).

    The magnetisation phase phi(x; temperature, no_of_sites) is defined via
    m(x; temperature, no_of_sites) = (|| m(x; temperature, no_of_sites) ||, phi(x; temperature, no_of_sites))^t in
    radial coordinates, where m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the
    Cartesian magnetisation, with x_i the position of particle i at the time of observation.

    NB, RuntimeWarnings are often thrown by lines 175, 178, 179 and 188 of scipy/stats/_hypotests.py when calculating
    the p-value of scipy.stats.cramervonmises() for highly nonergodic magnetisation-phase samples.  They do not affect
    the calculation of the Cramér-von Mises integral itself.

    Parameters
    ----------
    sample_directory : str
        The location of the directory containing the sample and Metropolis acceptance rate.
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.

    Returns
    -------
    numpy.ndarray
        The normalised Cramér-von Mises integral.  A float estimating the normalised Cramér-von Mises integral.
    """
    return stats.cramervonmises(sample_getter.get_magnetisation_phase(
        sample_directory, temperature, temperature_index, no_of_sites)[no_of_equilibration_sweeps:], cdf="uniform",
                                args=(-math.pi, 2.0 * math.pi)).statistic
