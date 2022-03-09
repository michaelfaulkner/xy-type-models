from scipy import stats
import importlib
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import time

# import additional modules; have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")
setup_scripts = importlib.import_module("setup_scripts")
sample_getter = importlib.import_module("sample_getter")
run_script = importlib.import_module("run")


def main():
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    linear_system_sizes = [2 ** (index + 3) for index in range(5)]
    config_files_local = [f"config_files/cramer_von_mises_figure/{value}x{value}_local_moves.txt"
                          for value in linear_system_sizes]
    config_file_all = "config_files/cramer_von_mises_figure/all_moves.txt"
    (algorithm_name, output_directory_local, _, no_of_equilibration_sweeps_local, no_of_observations_local,
     temperatures_local, use_external_global_moves_local, external_global_moves_string_local, no_of_jobs_local,
     max_no_of_cpus) = run_script.get_config_data(config_files_local[0])
    (_, output_directory_all, _, no_of_equilibration_sweeps_all, no_of_observations_all, temperatures_all,
     use_external_global_moves_local, external_global_moves_string_local, no_of_jobs_all,
     _) = run_script.get_config_data(config_file_all)
    try:
        sample_getter.get_acceptance_rates(f"{output_directory_local}/job_0", 0)
    except OSError:
        raise Exception("Local-dynamics simulations have not been run - enter 'python run.py "
                        "config_files/cramer_von_mises_figure/local_moves.txt' in the top directory.")
    '''try:
        sample_getter.get_acceptance_rates(f"{output_directory_all}/job_0", 0)
    except OSError:
        raise Exception("All-dynamics simulations have not been run - enter 'python run.py "
                        "config_files/cramer_von_mises_figure/local_moves.txt' in the top directory.")'''
    output_directory = output_directory_local.replace("/8x8_local_moves", "")
    pool = setup_scripts.setup_pool(no_of_jobs_local, max_no_of_cpus)

    figure, axis = plt.subplots(1)
    axis.tick_params(which='both', width=2)
    axis.tick_params(which='major', length=7, labelsize=18, pad=10)
    axis.tick_params(which='minor', length=4)
    [axis.spines[spine].set_linewidth(2) for spine in ["top", "bottom", "left", "right"]]
    axis.set_xlabel(r"$1 / (\beta J)$", fontsize=20, labelpad=8)
    axis.set_ylabel(r"$n \omega_n^2$", fontsize=20, labelpad=8)

    colors = ["black", "red", "blue", "green", "yellow"]

    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        try:
            with open(f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}_"
                      f"{external_global_moves_string_local}_{length}x{length}_sites_{no_of_observations_local}_"
                      f"obs.tsv", "r") as output_file:
                output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                    if not line.startswith('#')]).transpose()
                cvms, cvm_errors = output_file_sans_header[1], output_file_sans_header[2]
        except (IOError, IndexError) as _:
            cvm_file = open(
                f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}_"
                f"{external_global_moves_string_local}_{length}x{length}_sites_{no_of_observations_local}_obs.tsv", "w")
            cvm_file.write("# temperature".ljust(30) + "n omega_n^2".ljust(30) + "n omega_n^2 error" + "\n")
            cvms, cvm_errors = [], []
            for temperature_index, temperature in enumerate(temperatures_local):
                print(f"Temperature = {temperature:.4f}")
                cvm_vs_job = pool.starmap(get_cramervonmises_of_magnetisation_phase, [
                    (f"{output_directory}/{length}x{length}_local_moves/job_{job_index}", temperature,
                     temperature_index, length ** 2, no_of_equilibration_sweeps_local) for job_index in
                    range(no_of_jobs_local)])
                cvm, cvm_error = np.mean(cvm_vs_job), (np.std(cvm_vs_job) / len(cvm_vs_job) ** 0.5)
                cvms.append(cvm)
                cvm_errors.append(cvm_error)
                cvm_file.write(f"{temperature:.14e}".ljust(30) + f"{cvm:.14e}".ljust(30) + f"{cvm_error:.14e}" + "\n")
            cvm_file.close()
        axis.semilogy(temperatures_local, cvms, marker=".", markersize=10, color=colors[system_size_index],
                      linestyle="None", label=fr"$N$ = {length}x{length}")
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()
    legend = axis.legend(loc="upper right", fontsize=10)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(1.5)
    figure.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}.pdf",
                   bbox_inches="tight")


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


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: no positional arguments allowed.")
    else:
        main()