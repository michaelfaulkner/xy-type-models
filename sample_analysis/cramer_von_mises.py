"""This module contains methods that compute the Cramér-von Mises statistics of the magnetisation phase."""
from sample_getter import get_magnetisation_phase
from scipy import stats
import math
import numpy as np


def get_cvm_mag_phase_vs_temperature(algorithm_name, output_directory, sample_directory, no_of_equilibration_sweeps,
                                     no_of_samples, temperatures, external_global_moves_string, no_of_runs, pool,
                                     length):
    try:
        with open(f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}"
                  f"_{external_global_moves_string}_{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_"
                  f"{temperatures[-1]:.4f}_{no_of_samples}_obs_{no_of_runs}_runs.tsv", "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            cvms, cvm_errors = output_file_sans_header[1], output_file_sans_header[2]
    except IOError:
        cvm_file = open(f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}_"
                        f"{external_global_moves_string}_{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_"
                        f"{temperatures[-1]:.4f}_{no_of_samples}_obs_{no_of_runs}_runs.tsv", "w")
        cvm_file.write("# temperature".ljust(30) + "n omega_n^2".ljust(30) + "n omega_n^2 error" + "\n")
        cvms, cvm_errors = [], []
        for temperature_index, temperature in enumerate(temperatures):
            print(f"Temperature = {temperature:.4f}")
            cvm_vs_run = pool.starmap(get_cramervonmises_of_magnetisation_phase, [
                (f"{sample_directory}/run_{run_index}", temperature, temperature_index, length ** 2,
                 no_of_equilibration_sweeps) for run_index in range(no_of_runs)])
            cvm, cvm_error = np.mean(cvm_vs_run), np.std(cvm_vs_run) / len(cvm_vs_run) ** 0.5
            cvms.append(cvm), cvm_errors.append(cvm_error)
            cvm_file.write(f"{temperature:.14e}".ljust(30) + f"{cvm:.14e}".ljust(30) + f"{cvm_error:.14e}" + "\n")
        cvm_file.close()
    return cvms, cvm_errors


def get_cramervonmises_of_magnetisation_phase(sample_directory, temperature, temperature_index, no_of_sites,
                                              no_of_equilibration_sweeps=None, thinning_level=None):
    r"""
    Returns (for the sample of the magnetisation phase) the normalised Cramér-von Mises statistic
    n \omega_n^2 := n \int (F_n(x) - F(x))^2 dF(x), where n is the sample size and F(x) is the CDF of the continuous
    uniform distribution Uniform(-pi, pi).

    The magnetisation phase phi(x; temperature, no_of_sites) is defined via
    m(x; temperature, no_of_sites) = (|| m(x; temperature, no_of_sites) ||, phi(x; temperature, no_of_sites))^t in
    radial coordinates, where m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the
    Cartesian magnetisation, with x_i the position of particle i at the time of sample.

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
    no_of_equilibration_sweeps : None or int, optional
        The total number of equilibration iterations of the Markov process.  If None, the entire sample is returned.
    thinning_level : None or int, optional
        The number of samples to be discarded between retained samples of the thinning process.  If None,
        all samples are retained.

    Returns
    -------
    numpy.ndarray
        The normalised Cramér-von Mises integral.  A float estimating the normalised Cramér-von Mises integral.
    """
    return stats.cramervonmises(get_magnetisation_phase(sample_directory, temperature, temperature_index, no_of_sites,
                                                        no_of_equilibration_sweeps, thinning_level),
                                cdf="uniform", args=(-math.pi, 2.0 * math.pi)).statistic
