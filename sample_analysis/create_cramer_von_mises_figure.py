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
    linear_system_sizes = [2 ** (index + 3) for index in range(3)]
    config_files_metrop_local = [f"config_files/cvm_figure/{value}x{value}_metrop_local_moves.txt" for value in
                                 linear_system_sizes]
    config_files_ecmc = [f"config_files/cvm_figure/{value}x{value}_ecmc.txt" for value in
                         linear_system_sizes]
    config_files_metrop_all = [f"config_files/cvm_figure/{value}x{value}_metrop_all_moves.txt" for value in
                               linear_system_sizes]

    (algorithm_name_metrop, output_directory, _, no_of_equilibration_sweeps_metrop_local,
     no_of_observations_metrop_local, temperatures_cmv, use_external_global_moves_cmv, external_global_moves_string_cmv,
     no_of_jobs_metrop_local, max_no_of_cpus) = run_script.get_config_data(config_files_metrop_local[0])
    (algorithm_name_ecmc, _, _, no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, _, _, _, no_of_jobs_ecmc, _
     ) = run_script.get_config_data(config_files_ecmc[0])
    (_, _, _, no_of_equilibration_sweeps_metrop_all, no_of_observations_metrop_all, temperatures_metrop_all,
     use_external_global_moves_metrop_all, external_global_moves_string_metrop_all, no_of_jobs_metrop_all,
     _) = run_script.get_config_data(config_files_metrop_all[0])
    output_directory = output_directory.replace("/8x8_metrop_local_moves", "")

    pool = setup_scripts.setup_pool(no_of_jobs_metrop_local, max_no_of_cpus)
    figure, axis = plt.subplots(1)
    axis.tick_params(which='both', width=3)
    axis.tick_params(which='major', length=7, labelsize=18, pad=10)
    axis.tick_params(which='minor', length=4)
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]
    axis.set_xlabel(r"$1 / (\beta J)$", fontsize=20, labelpad=8)
    axis.set_ylabel(r"$n \omega_n^2$", fontsize=20, labelpad=8)
    axis.set_yscale('log')

    inset_axis = plt.axes([0.685, 0.66, 0.2, 0.2])
    inset_axis.set_xlabel(r"$1 / (\beta J)$", fontsize=8, labelpad=2)
    inset_axis.set_ylabel(r"$p(\rm{twist})$", fontsize=8, labelpad=2)
    inset_axis.tick_params(which='major', width=3, labelsize=8)
    [inset_axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]

    colors = ["black", "red", "blue", "green", "yellow", "cyan"]

    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        cvm_metrops, cvm_metrop_errors = get_cramer_von_mises_vs_temperature(
            algorithm_name_metrop, output_directory, f"{output_directory}/{length}x{length}_metrop_local_moves",
            no_of_equilibration_sweeps_metrop_local, no_of_observations_metrop_local, temperatures_cmv,
            external_global_moves_string_cmv, no_of_jobs_metrop_local, pool, length)
        cvm_ecmcs, cvm_ecmc_errors = get_cramer_von_mises_vs_temperature(
            algorithm_name_ecmc, output_directory, f"{output_directory}/{length}x{length}_ecmc",
            no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures_cmv,
            external_global_moves_string_cmv, no_of_jobs_ecmc, pool, length)

        try:
            with open(f"{output_directory}/probability_of_global_twists_{algorithm_name_metrop.replace('-', '_')}_"
                      f"{length}x{length}_sites.tsv", "r") as output_file:
                output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                    if not line.startswith('#')]).transpose()
                twist_probabilities, twist_probability_errors = output_file_sans_header[1], output_file_sans_header[2]
        except (IOError, IndexError) as _:
            twists_file = open(f"{output_directory}/probability_of_global_twists_"
                               f"{algorithm_name_metrop.replace('-', '_')}_{length}x{length}_sites.tsv", "w")
            twists_file.write("# temperature".ljust(30) + "p(twist)".ljust(30) + "p(twist) error" + "\n")
            twist_probabilities, twist_probability_errors = [], []
            for temperature_index, temperature in enumerate(temperatures_metrop_all):
                twist_probability_vs_job = pool.starmap(sample_getter.get_acceptance_rates, [
                    (f"{output_directory}/{length}x{length}_metrop_all_moves/job_{job_index}", temperature_index)
                    for job_index in range(no_of_jobs_metrop_all)])
                twist_probability, twist_probability_error = (np.mean(twist_probability_vs_job), np.std(
                    twist_probability_vs_job) / len(twist_probability_vs_job) ** 0.5)
                twist_probabilities.append(twist_probability)
                twist_probability_errors.append(twist_probability_error)
                twists_file.write(f"{temperature:.14e}".ljust(30) + f"{twist_probability:.14e}".ljust(30) +
                                  f"{twist_probability_error:.14e}" + "\n")
            twists_file.close()

        try:
            with open(f"{output_directory}/physical_time_steps_{algorithm_name_metrop.replace('-', '_')}_"
                      f"{external_global_moves_string_cmv}_{length}x{length}_sites.tsv", "r") as _:
                pass
        except IOError:
            physical_time_step_file = open(
                f"{output_directory}/physical_time_steps_{algorithm_name_metrop.replace('-', '_')}_"
                f"{external_global_moves_string_cmv}_{length}x{length}_sites.tsv", "w")
            physical_time_step_file.write("# temperature".ljust(30) + "Delta t".ljust(30) + "Delta t error" + "\n")
            for temperature_index, temperature in enumerate(temperatures_cmv):
                print(temperature_index)
                physical_time_step_vs_job = pool.starmap(sample_getter.get_physical_time_step, [
                    (algorithm_name_metrop, f"{output_directory}/{length}x{length}_metrop_local_moves/job_{job_index}",
                     temperature_index) for job_index in range(no_of_jobs_metrop_all)])
                physical_time_step_file.write(
                    f"{temperature:.14e}".ljust(30) + f"{np.mean(physical_time_step_vs_job):.14e}".ljust(30) +
                    f"{np.std(physical_time_step_vs_job) / len(physical_time_step_vs_job) ** 0.5:.14e}" + "\n")
            physical_time_step_file.close()

        axis.errorbar(temperatures_cmv, cvm_metrops, cvm_metrop_errors, marker=".", markersize=10,
                      color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length} (Metropolis)")
        axis.errorbar(temperatures_cmv, cvm_ecmcs, cvm_ecmc_errors, marker="*", markersize=10,
                      color=colors[system_size_index], linestyle="None",
                      label=fr"$N$ = {length}x{length} (event-chain)")
        inset_axis.errorbar(temperatures_metrop_all, twist_probabilities, twist_probability_errors, marker=".",
                            markersize=10, color=colors[system_size_index], linestyle="None")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()
    legend = axis.legend(loc="center left", fontsize=8)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(3)
    figure.tight_layout()
    figure.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_and_ecmc.pdf",
                   bbox_inches="tight")


def get_cramer_von_mises_vs_temperature(algorithm_name, output_directory, sample_directory, no_of_equilibration_sweeps,
                                        no_of_observations, temperatures, external_global_moves_string, no_of_jobs,
                                        pool, length):
    try:
        with open(f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}"
                  f"_{external_global_moves_string}_{length}x{length}_sites_{no_of_observations}_obs.tsv",
                  "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            cvms, cvm_errors = output_file_sans_header[1], output_file_sans_header[2]
    except (IOError, IndexError) as _:
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


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: no positional arguments allowed.")
    else:
        main()
