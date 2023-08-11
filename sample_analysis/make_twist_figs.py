from matplotlib.ticker import MultipleLocator
from sample_getter import get_acceptance_rates, get_magnetisation_phase, get_physical_time_step
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
    """n.b., the following mid-temps, lower-trans, upper-trans and high-temps sims are part of the CvM figure; we also 
        use them here for the twist-prob figure"""
    base_config_file_mid_temps = f"config_files/cvm_figs/4x4_metrop_low_temps.txt"
    base_config_file_lower_trans = f"config_files/cvm_figs/4x4_metrop_lower_trans.txt"
    base_config_file_upper_trans = f"config_files/cvm_figs/4x4_metrop_upper_trans.txt"
    base_config_file_high_temps = f"config_files/cvm_figs/4x4_metrop_high_temps.txt"

    (algorithm_name, sample_directory_4x4_low_temp_all, _, _, no_of_equilibration_sweeps_low_temp,
     no_of_observations_low_temp, temperatures_low_temp, _, external_global_moves_string_all, no_of_jobs_low_temp, _,
     max_no_of_cpus) = run_script.get_config_data(base_config_file_low_temp_all)
    external_global_moves_string_local = run_script.get_config_data(base_config_file_low_temp_local)[8]
    (_, sample_directory_4x4_mid_temps, _, _, no_of_equilibration_sweeps_mid_temps,
     no_of_observations_small_systems_mid_temps, temperatures_mid_temps, _, _, no_of_jobs_mid_temps, _, _
     ) = run_script.get_config_data(base_config_file_mid_temps)
    (_, _, _, _, no_of_equilibration_sweeps_lower_trans, no_of_observations_lower_trans, temperatures_lower_trans, _, _,
     no_of_jobs_lower_trans, _, _) = run_script.get_config_data(base_config_file_lower_trans)
    (_, _, _, _, no_of_equilibration_sweeps_upper_trans, no_of_observations_upper_trans, temperatures_upper_trans, _, _,
     no_of_jobs_upper_trans, _, _) = run_script.get_config_data(base_config_file_upper_trans)
    (_, _, _, _, no_of_equilibration_sweeps_high_temps, no_of_observations_high_temps, temperatures_high_temps, _, _,
     no_of_jobs_high_temps, _, _) = run_script.get_config_data(base_config_file_high_temps)
    temperatures = [*temperatures_low_temp, *temperatures_mid_temps, *temperatures_lower_trans,
                    *temperatures_upper_trans, *temperatures_high_temps]
    approx_transition_temperature = 0.887
    reduced_temperatures = [temperature / approx_transition_temperature for temperature in temperatures]

    """For L = 128, we split the four low-temperature CvM simulations (mid_temps for twist analysis) into two 
        two-temperature directories.  This kept the CvM simulations within the two-week time limit on BlueCrystal 4, as 
        this larger system size uses more CPU time at fixed simulation timescale, but also requires longer simulation 
        timescales for CvM convergence.  Due to these longer required simulation timescales, L = 64 and 128 both 
        require the additional field no_of_observations_large_systems_mid_temps."""
    base_config_file_128x128_mid_temps_lower = f"config_files/cvm_figs/128x128_metrop_lowest_temps.txt"
    base_config_file_128x128_mid_temps_upper = f"config_files/cvm_figs/128x128_metrop_low_temps.txt"
    temperatures_128x128_mid_temps_lower = run_script.get_config_data(base_config_file_128x128_mid_temps_lower)[6]
    (_, _, _, _, _, no_of_observations_large_systems_mid_temps, temperatures_128x128_mid_temps_upper, _, _, _,
     _, _) = run_script.get_config_data(base_config_file_128x128_mid_temps_upper)
    no_of_observations_mid_temps = [no_of_observations_small_systems_mid_temps if index < 4 else
                                    no_of_observations_large_systems_mid_temps for index in
                                    range(len(linear_system_sizes))]

    output_directory_low_temp = sample_directory_4x4_low_temp_all.replace("/4x4_low_temp_all", "")
    sample_directories_low_temp_all = [f"{output_directory_low_temp}/{length}x{length}_low_temp_all" for length in
                                       linear_system_sizes]
    sample_directories_low_temp_local = [f"{output_directory_low_temp}/{length}x{length}_low_temp_local" for length in
                                         linear_system_sizes]
    output_directory_higher_temps = sample_directory_4x4_mid_temps.replace("/4x4_metrop_low_temps", "")
    sample_directories_mid_temps = [f"{output_directory_higher_temps}/{length}x{length}_metrop_low_temps" for length in
                                    linear_system_sizes]
    sample_directories_lower_trans = [f"{output_directory_higher_temps}/{length}x{length}_metrop_lower_trans"
                                      for length in linear_system_sizes]
    sample_directories_upper_trans = [f"{output_directory_higher_temps}/{length}x{length}_metrop_upper_trans"
                                      for length in linear_system_sizes]
    sample_directories_high_temps = [f"{output_directory_higher_temps}/{length}x{length}_metrop_high_temps"
                                     for length in linear_system_sizes]
    """The following line defines an additional sample directory for L = 128, as the L = 128 CvM simulations 
        were split into two two-temperature directories.  This kept CvM simulations within the two-week time limit on 
        BlueCrystal 4, as this larger system size uses more CPU time at fixed simulation timescale, but also requires 
        longer simulation timescales for CvM convergence."""
    sample_directory_128x128_mid_temps_lower = f"{output_directory_higher_temps}/128x128_metrop_lowest_temps"
    pool = setup_pool(no_of_jobs_low_temp, max_no_of_cpus)

    fig, axes = plt.subplots(2, 1, figsize=(5.0 * 5.75 / 5, 4.5 * 5 / 5.75))
    fig.text(0.115, 0.875, "b", fontsize=20, weight='bold')
    fig.text(0.115, 0.3925, "c", fontsize=20, weight='bold')
    fig.tight_layout(h_pad=1.75)

    axes[0].set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=20, labelpad=-12)
    axes[0].set_ylabel(r"$p_{\rm twist}$", fontsize=20, labelpad=1)
    axes[0].set_yscale('log')
    axes[0].set_ylim([8.0 * 10 ** (-8), 1.5])
    axes[0].set_yticks([10 ** (-7), 10 ** (-5), 10 ** (-3), 10 ** (-1)])
    axes[1].set_xlabel(r"$1 / \ln N$", fontsize=20, labelpad=-0)
    axes[1].set_ylabel(r"$\frac{\langle s_{\phi_{\mathbf{m}}}^2 \left(\beta J \! = \! 10, \! n \! = \! 10^6 \right) "
                       r"\rangle}{{\rm Var}\left[\phi_{\mathbf{m}} \right]}$", fontsize=20, labelpad=6)
    axes[1].set_ylim([-0.1, 1.1])
    axes[1].yaxis.set_minor_locator(MultipleLocator(base=0.5))

    [axis.spines[spine].set_linewidth(3.0) for spine in ["top", "bottom", "left", "right"] for axis in axes]
    for axis in axes:
        axis.tick_params(which='major', direction='in', width=2.5, length=5, labelsize=18, pad=2.5)
        axis.tick_params(which='minor', direction='in', width=1.5, length=4)

    colors = ["black", "red", "blue", "green", "magenta", "indigo"][:no_of_system_sizes]
    colors.reverse()

    low_temp_simulation_variance_vs_system_size_local = []
    low_temp_simulation_variance_error_vs_system_size_local = []
    low_temp_simulation_variance_vs_system_size_all = []
    low_temp_simulation_variance_error_vs_system_size_all = []
    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        print(f"Number of sites, N = {length}x{length}")
        """compute non-used physical time steps for our records"""
        compute_physical_time_steps(algorithm_name, external_global_moves_string_all, output_directory_low_temp,
                                    sample_directories_low_temp_all[system_size_index], temperatures_low_temp, length,
                                    no_of_observations_low_temp, no_of_jobs_low_temp, pool)

        twist_probabilities_low_temp, twist_probability_errors_low_temp = get_twist_probabilities_and_errors(
            algorithm_name, output_directory_low_temp, sample_directories_low_temp_all[system_size_index],
            temperatures_low_temp, length, no_of_observations_low_temp, no_of_jobs_low_temp, pool)
        if length < 128:
            twist_probabilities_mid_temps, twist_probability_errors_mid_temps = get_twist_probabilities_and_errors(
                algorithm_name, output_directory_higher_temps, sample_directories_mid_temps[system_size_index],
                temperatures_mid_temps, length, no_of_observations_mid_temps[system_size_index], no_of_jobs_mid_temps,
                pool)
        else:
            twist_probs_mid_temps_lower, twist_prob_errors_mid_temps_lower = get_twist_probabilities_and_errors(
                algorithm_name, output_directory_higher_temps, sample_directory_128x128_mid_temps_lower,
                temperatures_128x128_mid_temps_lower, length, no_of_observations_mid_temps[system_size_index],
                no_of_jobs_mid_temps, pool)
            twist_probs_mid_temps_upper, twist_prob_errors_mid_temps_upper = get_twist_probabilities_and_errors(
                algorithm_name, output_directory_higher_temps, sample_directories_mid_temps[system_size_index],
                temperatures_128x128_mid_temps_upper, length, no_of_observations_mid_temps[system_size_index],
                no_of_jobs_mid_temps, pool)
            twist_probabilities_mid_temps = [*twist_probs_mid_temps_lower, *twist_probs_mid_temps_upper]
            twist_probability_errors_mid_temps = [*twist_prob_errors_mid_temps_lower,
                                                  *twist_prob_errors_mid_temps_upper]
        twist_probabilities_lower_trans, twist_probability_errors_lower_trans = get_twist_probabilities_and_errors(
            algorithm_name, output_directory_higher_temps, sample_directories_lower_trans[system_size_index],
            temperatures_lower_trans, length, no_of_observations_lower_trans, no_of_jobs_lower_trans, pool)
        twist_probabilities_upper_trans, twist_probability_errors_upper_trans = get_twist_probabilities_and_errors(
            algorithm_name, output_directory_higher_temps, sample_directories_upper_trans[system_size_index],
            temperatures_upper_trans, length, no_of_observations_upper_trans, no_of_jobs_upper_trans, pool)
        twist_probabilities_high_temps, twist_probability_errors_high_temps = get_twist_probabilities_and_errors(
            algorithm_name, output_directory_higher_temps, sample_directories_high_temps[system_size_index],
            temperatures_high_temps, length, no_of_observations_high_temps, no_of_jobs_high_temps, pool)
        twist_probabilities = [*twist_probabilities_low_temp, *twist_probabilities_mid_temps,
                               *twist_probabilities_lower_trans, *twist_probabilities_upper_trans,
                               *twist_probabilities_high_temps]
        twist_probability_errors = [*twist_probability_errors_low_temp, *twist_probability_errors_mid_temps,
                                    *twist_probability_errors_lower_trans, *twist_probability_errors_upper_trans,
                                    *twist_probability_errors_high_temps]

        """n.b., the 1 / (math.pi ** 2 / 3.0) factors below divide the simulation variances by their expected value"""
        (low_temp_simulation_variances_local,
         low_temp_simulation_variance_errors_local) = np.array(get_mag_phase_simulation_variances_and_errors(
            algorithm_name, external_global_moves_string_local, output_directory_low_temp,
            sample_directories_low_temp_local[system_size_index], temperatures_low_temp, length,
            no_of_equilibration_sweeps_low_temp, no_of_observations_low_temp, no_of_jobs_low_temp, pool)) / (
                math.pi ** 2 / 3.0)
        (low_temp_simulation_variances_all,
         low_temp_simulation_variance_errors_all) = np.array(get_mag_phase_simulation_variances_and_errors(
            algorithm_name, external_global_moves_string_all, output_directory_low_temp,
            sample_directories_low_temp_all[system_size_index], temperatures_low_temp, length,
            no_of_equilibration_sweeps_low_temp, no_of_observations_low_temp, no_of_jobs_low_temp, pool)) / (
                math.pi ** 2 / 3.0)

        low_temp_simulation_variance_vs_system_size_local.append(low_temp_simulation_variances_local[0])
        low_temp_simulation_variance_error_vs_system_size_local.append(low_temp_simulation_variance_errors_local[0])
        low_temp_simulation_variance_vs_system_size_all.append(low_temp_simulation_variances_all[0])
        low_temp_simulation_variance_error_vs_system_size_all.append(low_temp_simulation_variance_errors_all[0])

        axes[0].errorbar(reduced_temperatures, twist_probabilities, twist_probability_errors, marker=".", markersize=11,
                         color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length}")

    inverse_log_system_sizes = [1.0 / np.log(length ** 2) for length in linear_system_sizes]
    axes[1].errorbar(inverse_log_system_sizes, low_temp_simulation_variance_vs_system_size_local,
                     low_temp_simulation_variance_error_vs_system_size_local, marker="*", markersize=11, color="black",
                     linestyle='None', label="local Metropolis dynamics")
    axes[1].errorbar(inverse_log_system_sizes, low_temp_simulation_variance_vs_system_size_all,
                     low_temp_simulation_variance_error_vs_system_size_all, marker=".", markersize=11, color="red",
                     linestyle='None', label="local Metropolis dynamics with twists")

    legends = [axes[0].legend(loc="lower right", ncol=2, fontsize=10.5),
               axes[1].legend(loc="lower right", fontsize=10.5)]
    [legend.get_frame().set_edgecolor("k") for legend in legends]
    [legend.get_frame().set_lw(3) for legend in legends]
    fig.savefig(f"{output_directory_low_temp}/twist_probability_vs_temperature_and_magnetisation_phase_simulation_"
                f"variance_vs_system_size_xy_gaussian_noise_metropolis.pdf", bbox_inches="tight")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()


def compute_physical_time_steps(algorithm_name, external_global_moves_string, output_directory, sample_directory,
                                temperatures, length, no_of_observations, no_of_jobs, pool):
    try:
        with open(f"{output_directory}/physical_time_steps_{algorithm_name.replace('-', '_')}_"
                  f"{external_global_moves_string}_{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_"
                  f"{temperatures[-1]:.4f}_{no_of_observations}_obs_{no_of_jobs}_jobs.tsv", "r") as _:
            pass
    except IOError:
        physical_time_step_file = open(
            f"{output_directory}/physical_time_steps_{algorithm_name.replace('-', '_')}_{external_global_moves_string}_"
            f"{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_{temperatures[-1]:.4f}_{no_of_observations}_"
            f"obs_{no_of_jobs}_jobs.tsv", "w")
        physical_time_step_file.write("# temperature".ljust(30) + "Delta t".ljust(30) + "Delta t error" + "\n")
        for temperature_index, temperature in enumerate(temperatures):
            physical_time_step_vs_job = pool.starmap(get_physical_time_step, [
                (algorithm_name, f"{sample_directory}/job_{job_index}", temperature_index) for job_index in
                range(no_of_jobs)])
            physical_time_step_file.write(
                f"{temperature:.14e}".ljust(30) + f"{np.mean(physical_time_step_vs_job):.14e}".ljust(30) +
                f"{np.std(physical_time_step_vs_job) / len(physical_time_step_vs_job) ** 0.5:.14e}" + "\n")
        physical_time_step_file.close()


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


def get_mag_phase_simulation_variances_and_errors(algorithm_name, external_global_moves_string, output_directory,
                                                  sample_directory, temperatures, length, no_of_equilibration_sweeps,
                                                  no_of_observations, no_of_jobs, pool):
    try:
        with open(f"{output_directory}/magnetisation_phase_simulation_variance_{algorithm_name.replace('-', '_')}_"
                  f"{external_global_moves_string}_{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_"
                  f"{temperatures[-1]:.4f}_{no_of_observations}_obs_{no_of_jobs}_jobs.tsv", "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            simulation_variances, simulation_variance_errors = output_file_sans_header[1], output_file_sans_header[2]
    except IOError:
        simulation_variance_file = open(
            f"{output_directory}/magnetisation_phase_simulation_variance_{algorithm_name.replace('-', '_')}_"
            f"{external_global_moves_string}_{length}x{length}_sites_temp_range_{temperatures[0]:.4f}_to_"
            f"{temperatures[-1]:.4f}_{no_of_observations}_obs_{no_of_jobs}_jobs.tsv", "w")
        simulation_variance_file.write("# temperature".ljust(30) + "simulation variance".ljust(30) +
                                       "simulation variance error" + "\n")
        simulation_variances, simulation_variance_errors = [], []
        for temperature_index, temperature in enumerate(temperatures):
            simulation_variance_vs_job = pool.starmap(get_simulation_variance_of_magnetisation_phase, [
                (f"{sample_directory}/job_{job_index}", temperature, temperature_index, length ** 2,
                 no_of_equilibration_sweeps) for job_index in range(no_of_jobs)])
            simulation_variance, simulation_variance_error = np.mean(simulation_variance_vs_job), np.std(
                simulation_variance_vs_job) / len(simulation_variance_vs_job) ** 0.5
            simulation_variances.append(simulation_variance)
            simulation_variance_errors.append(simulation_variance_error)
            simulation_variance_file.write(f"{temperature:.14e}".ljust(30) + f"{simulation_variance:.14e}".ljust(30) +
                                           f"{simulation_variance_error:.14e}" + "\n")
        simulation_variance_file.close()
    return simulation_variances, simulation_variance_errors


def get_simulation_variance_of_magnetisation_phase(sample_directory, temperature, temperature_index, no_of_sites,
                                                   no_of_equilibration_sweeps):
    r"""
    Returns the simulation variance of the magnetisation phase.

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
        The simulation variance of the magnetisation phase.  A float.
    """
    return np.var(get_magnetisation_phase(sample_directory, temperature, temperature_index, no_of_sites)[
                  no_of_equilibration_sweeps:], ddof=1)


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("InterfaceError: At most one positional arguments permitted.  None are required but you may "
                        "provide no_of_system_sizes, which must be an integer greater than 0 and less than 7 (default "
                        "value is 6).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be no_of_system_sizes - which must be an integer greater "
              "than 0 and less than 7 (default value is 6).")
        chosen_no_of_system_sizes = int(sys.argv[1])
        if chosen_no_of_system_sizes < 1 or chosen_no_of_system_sizes > 7:
            raise Exception("InterfaceError: no_of_system_sizes must be an integer greater than 0 and less than 7 "
                            "(default value is 6).")
        main(chosen_no_of_system_sizes)
    else:
        print("No positional arguments provided.  None are required but you may provide no_of_system_sizes, which must "
              "be an integer greater than 0 and less than 7 (default value is 6).")
        main()
