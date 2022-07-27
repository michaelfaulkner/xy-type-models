from make_cvm_figs import compute_physical_time_steps
from sample_getter import get_acceptance_rates, get_magnetisation_phase
from setup_scripts import setup_pool
import importlib
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import time

# import run script - have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def main(no_of_system_sizes=6):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    linear_system_sizes = [2 ** (index + 2) for index in range(no_of_system_sizes)]
    base_config_file_low_temp_all = f"config_files/twist_figs/4x4_low_temp_all.txt"
    base_config_file_low_temp_local = f"config_files/twist_figs/4x4_low_temp_local.txt"
    """n.b., the following mid- and high-temps sims are part of the CvM figure; we also use them here for the 
        twist-prob figure"""
    base_config_file_mid_temps = f"config_files/cvm_fig/4x4_metrop_low_temps.txt"
    base_config_file_high_temps = f"config_files/cvm_fig/4x4_metrop_high_temps.txt"

    (algorithm_name, sample_directory_4x4_low_temp_all, _, _, no_of_equilibration_sweeps_low_temp,
     no_of_observations_low_temp, temperatures_low_temp, _, external_global_moves_string_all, no_of_jobs_low_temp, _,
     max_no_of_cpus) = run_script.get_config_data(base_config_file_low_temp_all)
    external_global_moves_string_local = run_script.get_config_data(base_config_file_low_temp_local)[8]
    (_, sample_directory_4x4_mid_temps, _, _, no_of_equilibration_sweeps_mid_temps,
     no_of_observations_mid_temps, temperatures_mid_temps, _, _, no_of_jobs_mid_temps, _, _
     ) = run_script.get_config_data(base_config_file_mid_temps)
    (_, _, _, _, no_of_equilibration_sweeps_high_temps, no_of_observations_high_temps, temperatures_high_temps, _, _,
     no_of_jobs_high_temps, _, _) = run_script.get_config_data(base_config_file_high_temps)
    temperatures = [*temperatures_low_temp, *temperatures_mid_temps, *temperatures_high_temps]
    approx_transition_temperature = 0.887
    reduced_temperatures = [temperature / approx_transition_temperature for temperature in temperatures]

    output_directory_low_temp = sample_directory_4x4_low_temp_all.replace("/4x4_low_temp_all", "")
    sample_directories_low_temp_all = [f"{output_directory_low_temp}/{length}x{length}_low_temp_all" for length in
                                       linear_system_sizes]
    sample_directories_low_temp_local = [f"{output_directory_low_temp}/{length}x{length}_low_temp_local" for length in
                                         linear_system_sizes]
    output_directory_higher_temps = sample_directory_4x4_mid_temps.replace("/4x4_metrop_low_temps", "")
    sample_directories_mid_temps = [f"{output_directory_higher_temps}/{length}x{length}_metrop_low_temps" for length in
                                    linear_system_sizes]
    sample_directories_high_temps = [f"{output_directory_higher_temps}/{length}x{length}_metrop_high_temps"
                                     for length in linear_system_sizes]
    pool = setup_pool(no_of_jobs_low_temp, max_no_of_cpus)

    fig, axes = plt.subplots(2, 1, figsize=(5.0, 4.5))
    fig.text(0.0025, 0.925, "f", fontsize=20, weight='bold')
    fig.text(0.0025, 0.4625, "g", fontsize=20, weight='bold')
    fig.tight_layout(h_pad=2.5)

    axes[0].set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=20, labelpad=-6)
    axes[0].set_yscale('log')
    axes[0].set_ylim([9.0 * 10 ** (-8), 0.7])
    axes[0].set_yticks([10 ** (-7), 10 ** (-5), 10 ** (-3), 10 ** (-1)])
    axes[0].set_ylabel(r"$p_{\rm twist}$", fontsize=20, labelpad=4)
    axes[1].set_xlabel(r"$N^{-1 / 2}$", fontsize=20, labelpad=3)
    axes[1].set_ylabel(r"$\frac{s_{\phi_m, n=10^6}^2(\beta J = 10)}{{\rm Var}\left[ \phi_m \right]}$", fontsize=20,
                       labelpad=1)
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes]
    for axis in axes:
        axis.tick_params(which='both', direction='in', width=3)
        axis.tick_params(which='major', length=5, labelsize=18, pad=5)
        axis.tick_params(which='minor', length=4)

    colors = ["black", "red", "blue", "green", "tab:brown", "magenta", "indigo"][:no_of_system_sizes]
    colors.reverse()

    low_temp_sample_variance_vs_system_size_local = []
    low_temp_sample_variance_error_vs_system_size_local = []
    low_temp_sample_variance_vs_system_size_all = []
    low_temp_sample_variance_error_vs_system_size_all = []
    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        print(f"Number of sites, N = {length}x{length}")
        compute_physical_time_steps(algorithm_name, external_global_moves_string_all, output_directory_low_temp,
                                    sample_directories_low_temp_all[system_size_index], temperatures_low_temp, length,
                                    no_of_observations_low_temp, no_of_jobs_low_temp, pool)

        twist_probabilities_low_temp, twist_probability_errors_low_temp = get_twist_probabilities_and_errors(
            algorithm_name, output_directory_low_temp, sample_directories_low_temp_all[system_size_index],
            temperatures_low_temp, length, no_of_observations_low_temp, no_of_jobs_low_temp, pool)
        twist_probabilities_mid_temps, twist_probability_errors_mid_temps = get_twist_probabilities_and_errors(
            algorithm_name, output_directory_low_temp, sample_directories_mid_temps[system_size_index],
            temperatures_mid_temps, length, no_of_observations_mid_temps, no_of_jobs_mid_temps, pool)
        twist_probabilities_high_temps, twist_probability_errors_high_temps = get_twist_probabilities_and_errors(
            algorithm_name, output_directory_low_temp, sample_directories_high_temps[system_size_index],
            temperatures_high_temps, length, no_of_observations_high_temps, no_of_jobs_high_temps, pool)
        twist_probabilities = [*twist_probabilities_low_temp, *twist_probabilities_mid_temps,
                               *twist_probabilities_high_temps]
        twist_probability_errors = [*twist_probability_errors_low_temp, *twist_probability_errors_mid_temps,
                                    *twist_probability_errors_high_temps]

        """n.b., the 1 / (math.pi ** 2 / 3.0) factors below divide the sample variances by their expected value"""
        (low_temp_sample_variances_local,
         low_temp_sample_variance_errors_local) = np.array(get_mag_phase_sample_variances_and_errors(
            algorithm_name, external_global_moves_string_local, output_directory_low_temp,
            sample_directories_low_temp_local[system_size_index], temperatures_low_temp, length,
            no_of_equilibration_sweeps_low_temp, no_of_observations_low_temp, no_of_jobs_low_temp, pool)) / (
                math.pi ** 2 / 3.0)
        (low_temp_sample_variances_all,
         low_temp_sample_variance_errors_all) = np.array(get_mag_phase_sample_variances_and_errors(
            algorithm_name, external_global_moves_string_all, output_directory_low_temp,
            sample_directories_low_temp_all[system_size_index], temperatures_low_temp, length,
            no_of_equilibration_sweeps_low_temp, no_of_observations_low_temp, no_of_jobs_low_temp, pool)) / (
                math.pi ** 2 / 3.0)
        _, _ = np.array(get_mag_phase_sample_variances_and_errors(
            algorithm_name, external_global_moves_string_all, output_directory_low_temp,
            sample_directories_mid_temps[system_size_index], temperatures_mid_temps, length,
            no_of_equilibration_sweeps_mid_temps, no_of_observations_mid_temps, no_of_jobs_mid_temps, pool)) / (
                math.pi ** 2 / 3.0)
        _, _ = np.array(get_mag_phase_sample_variances_and_errors(
            algorithm_name, external_global_moves_string_all, output_directory_low_temp,
            sample_directories_high_temps[system_size_index], temperatures_high_temps, length,
            no_of_equilibration_sweeps_high_temps, no_of_observations_high_temps, no_of_jobs_high_temps, pool)) / (
                math.pi ** 2 / 3.0)

        low_temp_sample_variance_vs_system_size_local.append(low_temp_sample_variances_local[0])
        low_temp_sample_variance_error_vs_system_size_local.append(low_temp_sample_variance_errors_local[0])
        low_temp_sample_variance_vs_system_size_all.append(low_temp_sample_variances_all[0])
        low_temp_sample_variance_error_vs_system_size_all.append(low_temp_sample_variance_errors_all[0])

        axes[0].errorbar(reduced_temperatures, twist_probabilities, twist_probability_errors, marker=".", markersize=10,
                         color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length}")

    inverse_linear_system_sizes = [1.0 / length for length in linear_system_sizes]
    axes[1].errorbar(inverse_linear_system_sizes, low_temp_sample_variance_vs_system_size_local,
                     low_temp_sample_variance_error_vs_system_size_local, marker=".", markersize=8, color="black",
                     linestyle='None', label="local dynamics only")
    axes[1].errorbar(inverse_linear_system_sizes, low_temp_sample_variance_vs_system_size_all,
                     low_temp_sample_variance_error_vs_system_size_all, marker="*", markersize=8, color="red",
                     linestyle='None', label="local dynamics with twists")

    legends = [axes[0].legend(loc="lower right", fontsize=10), axes[1].legend(loc="lower right", fontsize=10)]
    [legend.get_frame().set_edgecolor("k") for legend in legends]
    [legend.get_frame().set_lw(3) for legend in legends]
    fig.savefig(f"{output_directory_low_temp}/twist_probability_vs_temperature_and_magnetisation_phase_sample_variance_"
                f"vs_system_size_xy_gaussian_noise_metropolis_and_ecmc.pdf", bbox_inches="tight")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()


def get_twist_probabilities_and_errors(algorithm_name, output_directory, sample_directory, temperatures, length,
                                       no_of_observations, no_of_jobs, pool):
    try:
        with open(f"{output_directory}/probability_of_global_twists_{algorithm_name.replace('-', '_')}_"
                  f"{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_{temperatures[-1]:.4f}_"
                  f"{no_of_observations}_obs_{no_of_jobs}_jobs.tsv", "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            twist_probabilities, twist_probability_errors = output_file_sans_header[1], output_file_sans_header[2]
    except IOError:
        twists_file = open(
            f"{output_directory}/probability_of_global_twists_{algorithm_name.replace('-', '_')}_{length}x{length}_"
            f"sites_temp_range_{temperatures[0]:.4f}_to_{temperatures[-1]:.4f}_{no_of_observations}_obs_{no_of_jobs}_"
            f"jobs.tsv", "w")
        twists_file.write("# temperature".ljust(30) + "p(twist)".ljust(30) + "p(twist) error" + "\n")
        twist_probabilities, twist_probability_errors = [], []
        for temperature_index, temperature in enumerate(temperatures):
            twist_probability_vs_job = np.array(pool.starmap(get_acceptance_rates, [
                (f"{sample_directory}/job_{job_index}", temperature_index) for job_index in range(no_of_jobs)]))[:, 2]
            twist_probability, twist_probability_error = (np.mean(twist_probability_vs_job), np.std(
                twist_probability_vs_job) / len(twist_probability_vs_job) ** 0.5)
            twist_probabilities.append(twist_probability)
            twist_probability_errors.append(twist_probability_error)
            twists_file.write(f"{temperature:.14e}".ljust(30) + f"{twist_probability:.14e}".ljust(30) +
                              f"{twist_probability_error:.14e}" + "\n")
        twists_file.close()
    return twist_probabilities, twist_probability_errors


def get_mag_phase_sample_variances_and_errors(algorithm_name, external_global_moves_string, output_directory,
                                              sample_directory, temperatures, length, no_of_equilibration_sweeps,
                                              no_of_observations, no_of_jobs, pool):
    try:
        with open(f"{output_directory}/magnetisation_phase_sample_variance_{algorithm_name.replace('-', '_')}_"
                  f"{external_global_moves_string}_{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_"
                  f"{temperatures[-1]:.4f}_{no_of_observations}_obs_{no_of_jobs}_jobs.tsv", "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            sample_variances, sample_variance_errors = output_file_sans_header[1], output_file_sans_header[2]
    except IOError:
        sample_variance_file = open(
            f"{output_directory}/magnetisation_phase_sample_variance_{algorithm_name.replace('-', '_')}_"
            f"{external_global_moves_string}_{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_"
            f"{temperatures[-1]:.4f}_{no_of_observations}_obs_{no_of_jobs}_jobs.tsv", "w")
        sample_variance_file.write("# temperature".ljust(30) + "sample variance".ljust(30) + "sample variance error" +
                                   "\n")
        sample_variances, sample_variance_errors = [], []
        for temperature_index, temperature in enumerate(temperatures):
            sample_variance_vs_job = pool.starmap(get_sample_variance_of_magnetisation_phase, [
                (f"{sample_directory}/job_{job_index}", temperature, 0, length ** 2, no_of_equilibration_sweeps)
                for job_index in range(no_of_jobs)])
            sample_variance, sample_variance_error = np.mean(sample_variance_vs_job), np.std(
                sample_variance_vs_job) / len(sample_variance_vs_job) ** 0.5
            sample_variances.append(sample_variance)
            sample_variance_errors.append(sample_variance_error)
            sample_variance_file.write(f"{temperature:.14e}".ljust(30) + f"{sample_variance:.14e}".ljust(30) +
                                       f"{sample_variance_error:.14e}" + "\n")
        sample_variance_file.close()
    return sample_variances, sample_variance_errors


def get_sample_variance_of_magnetisation_phase(sample_directory, temperature, temperature_index, no_of_sites,
                                               no_of_equilibration_sweeps):
    r"""
    Returns the sample variance of the magnetisation phase.

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
        The sample variance of the magnetisation phase.  A float.
    """
    return np.var(get_magnetisation_phase(sample_directory, temperature, temperature_index, no_of_sites)[
                  no_of_equilibration_sweeps:])


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("InterfaceError: At most one positional arguments permitted.  None are required but you may "
                        "provide no_of_system_sizes, which must be an integer greater than 0 and less than 8 (default "
                        "value is 6).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be no_of_system_sizes - which must be an integer greater "
              "than 0 and less than 8 (default value is 6).")
        chosen_no_of_system_sizes = int(sys.argv[1])
        if chosen_no_of_system_sizes < 1 or chosen_no_of_system_sizes > 7:
            raise Exception("InterfaceError: no_of_system_sizes must be an integer greater than 0 and less than 8 "
                            "(default value is 6).")
        main(chosen_no_of_system_sizes)
    else:
        print("No positional arguments provided.  None are required but you may provide no_of_system_sizes, which must "
              "be an integer greater than 0 and less than 8 (default value is 6).")
        main()
