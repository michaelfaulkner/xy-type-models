from cramer_von_mises import get_cvm_mag_phase_vs_temperature
from sample_getter import get_physical_time_step
from setup_scripts import setup_pool
import importlib
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
    base_config_file_metrop_low_temps = f"config_files/cvm_figs/4x4_metrop_low_temps.txt"
    base_config_file_metrop_high_temps = f"config_files/cvm_figs/4x4_metrop_high_temps.txt"
    base_config_file_ecmc = f"config_files/cvm_figs/4x4_ecmc.txt"

    """temperatures_ecmc = [*temperatures_metrop_low_temps, *temperatures_metrop_high_temps] in the following lines"""
    (algorithm_name_metrop, sample_directory_4x4_metrop_low_temps, _, _, no_of_equilibration_sweeps_metrop_low_temps,
     no_of_observations_metrop_low_temps, temperatures_metrop_low_temps, _, external_global_moves_string,
     no_of_jobs_metrop_low_temps, _, max_no_of_cpus) = run_script.get_config_data(base_config_file_metrop_low_temps)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_high_temps, no_of_observations_metrop_high_temps,
     temperatures_metrop_high_temps, _, _, no_of_jobs_metrop_high_temps, _, _) = run_script.get_config_data(
        base_config_file_metrop_high_temps)
    (algorithm_name_ecmc, _, _, _, no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures_ecmc, _, _,
     no_of_jobs_ecmc, _, _) = run_script.get_config_data(base_config_file_ecmc)

    output_directory = sample_directory_4x4_metrop_low_temps.replace("/4x4_metrop_low_temps", "")
    sample_directories_metrop_low_temps = [f"{output_directory}/{length}x{length}_metrop_low_temps" for length in
                                           linear_system_sizes]
    sample_directories_metrop_high_temps = [f"{output_directory}/{length}x{length}_metrop_high_temps" for length in
                                            linear_system_sizes]
    sample_directories_ecmc = [f"{output_directory}/{length}x{length}_ecmc" for length in linear_system_sizes]
    pool = setup_pool(no_of_jobs_metrop_low_temps, max_no_of_cpus)

    approx_transition_temperature = 0.887
    reduced_temperatures_metrop_low_temps = [temperature / approx_transition_temperature for temperature in
                                             temperatures_metrop_low_temps]
    reduced_temperatures_metrop_high_temps = [temperature / approx_transition_temperature for temperature in
                                              temperatures_metrop_high_temps]
    reduced_temperatures_ecmc = [temperature / approx_transition_temperature for temperature in temperatures_ecmc]

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), gridspec_kw={'width_ratios': [1.2, 1.0]})
    fig.text(0.0025, 0.925, "d", fontsize=20, weight='bold')
    fig.text(0.52, 0.925, "e", fontsize=20, weight='bold')
    fig.tight_layout(w_pad=3.0)
    axes[0].set_yscale('log')
    axes[0].set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=20, labelpad=-0.5)
    axes[0].set_ylabel(r"$n \omega_{\phi_m,n}^2$", fontsize=20, labelpad=-7.5)
    axes[1].set_xlabel(r"$N^{-1 / 2}$", fontsize=20, labelpad=4)
    axes[1].set_ylabel(r"$\widetilde{\beta}_{\rm BKT} / \beta_{\rm int}$", fontsize=20, labelpad=-0.5)
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes]
    for axis_index, axis in enumerate(axes):
        axis.tick_params(which='both', direction='in', width=3)
        axis.tick_params(which='major', length=5, labelsize=18, pad=5)
        axis.tick_params(which='minor', length=4)

    inset_axis = plt.axes([0.6, 0.65, 0.16, 0.275])
    inset_axis.tick_params(which='both', direction='in', length=4, width=2, labelsize=12)
    inset_axis.yaxis.set_label_position("right"), inset_axis.yaxis.tick_right()
    inset_axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=12, labelpad=-0.5)
    inset_axis.set_ylabel(r"$n \omega_{\phi_m,n}^2$", fontsize=12, labelpad=-0.5)
    inset_axis.set_xlim([0.86, 1.325]), inset_axis.set_ylim([2.0, 5.0 * 10 ** 3])
    inset_axis.set_yscale('log')
    [inset_axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]

    colors = ["black", "red", "blue", "green", "tab:brown", "magenta", "indigo"][:no_of_system_sizes]
    colors.reverse()

    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        print(f"Number of sites, N = {length}x{length}")
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_low_temps[system_size_index], temperatures_metrop_low_temps, length,
            no_of_observations_metrop_low_temps, no_of_jobs_metrop_low_temps, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_high_temps[system_size_index], temperatures_metrop_high_temps, length,
            no_of_observations_metrop_high_temps, no_of_jobs_metrop_high_temps, pool)

        cvms_metrop_low_temps, cvm_errors_metrop_low_temps = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_low_temps[system_size_index],
            no_of_equilibration_sweeps_metrop_low_temps, no_of_observations_metrop_low_temps,
            temperatures_metrop_low_temps, external_global_moves_string, no_of_jobs_metrop_low_temps, pool, length)
        cvms_metrop_high_temps, cvm_errors_metrop_high_temps = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_high_temps[system_size_index],
            no_of_equilibration_sweeps_metrop_high_temps, no_of_observations_metrop_high_temps,
            temperatures_metrop_high_temps, external_global_moves_string, no_of_jobs_metrop_high_temps, pool, length)
        cvms_ecmc, cvm_errors_ecmc = get_cvm_mag_phase_vs_temperature(
            algorithm_name_ecmc, output_directory, sample_directories_ecmc[system_size_index],
            no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures_ecmc,
            external_global_moves_string, no_of_jobs_ecmc, pool, length)

        axes[0].errorbar(reduced_temperatures_metrop_low_temps, cvms_metrop_low_temps, cvm_errors_metrop_low_temps,
                         marker=".", markersize=10, color=colors[system_size_index], linestyle="None",
                         label=fr"$N$ = {length}x{length} Metrop.")
        axes[0].errorbar(reduced_temperatures_metrop_high_temps, cvms_metrop_high_temps, cvm_errors_metrop_high_temps,
                         marker=".", markersize=10, color=colors[system_size_index], linestyle="None")
        axes[0].errorbar(reduced_temperatures_ecmc, cvms_ecmc, cvm_errors_ecmc, marker="*", markersize=8,
                         color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length} ECMC")
        inset_axis.errorbar(reduced_temperatures_metrop_low_temps, cvms_metrop_low_temps, cvm_errors_metrop_low_temps,
                            marker=".", markersize=8, color=colors[system_size_index], linestyle="None")
        inset_axis.errorbar(reduced_temperatures_metrop_high_temps, cvms_metrop_high_temps, cvm_errors_metrop_high_temps,
                            marker=".", markersize=8, color=colors[system_size_index], linestyle="None")

    inverse_linear_system_sizes = [1.0 / length for length in linear_system_sizes]

    legends = [axes[0].legend(loc="upper right", fontsize=9)]  #, axes[1].legend(loc="upper left", fontsize=10)]
    [legend.get_frame().set_edgecolor("k") for legend in legends]
    [legend.get_frame().set_lw(3) for legend in legends]
    fig.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_and_ecmc.pdf",
                bbox_inches="tight")

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
