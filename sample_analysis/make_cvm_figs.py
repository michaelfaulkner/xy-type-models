from cramer_von_mises import get_cvm_mag_phase_vs_temperature
from make_twist_figs import compute_physical_time_steps, get_twist_probabilities_and_errors
from setup_scripts import setup_pool
import importlib
import matplotlib
import matplotlib.pyplot as plt
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
    base_config_file_metrop_lower_trans = f"config_files/cvm_figs/4x4_metrop_lower_trans.txt"
    base_config_file_metrop_upper_trans = f"config_files/cvm_figs/4x4_metrop_upper_trans.txt"
    base_config_file_metrop_high_temps = f"config_files/cvm_figs/4x4_metrop_high_temps.txt"
    base_config_file_ecmc = f"config_files/cvm_figs/4x4_ecmc.txt"

    """temperatures_ecmc = [*temperatures_metrop_low_temps, *temperatures_metrop_high_temps] in the following lines"""
    (algorithm_name_metrop, sample_directory_4x4_metrop_low_temps, _, _, no_of_equilibration_sweeps_metrop_low_temps,
     no_of_observations_metrop_low_temps, temperatures_metrop_low_temps, _, external_global_moves_string,
     no_of_jobs_metrop_low_temps, _, max_no_of_cpus) = run_script.get_config_data(base_config_file_metrop_low_temps)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_lower_trans, no_of_observations_metrop_lower_trans,
     temperatures_metrop_lower_trans, _, _, no_of_jobs_metrop_lower_trans, _, _) = run_script.get_config_data(
        base_config_file_metrop_lower_trans)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_upper_trans, no_of_observations_metrop_upper_trans,
     temperatures_metrop_upper_trans, _, _, no_of_jobs_metrop_upper_trans, _, _) = run_script.get_config_data(
        base_config_file_metrop_upper_trans)
    (_, _, _, _, no_of_equilibration_sweeps_metrop_high_temps, no_of_observations_metrop_high_temps,
     temperatures_metrop_high_temps, _, _, no_of_jobs_metrop_high_temps, _, _) = run_script.get_config_data(
        base_config_file_metrop_high_temps)
    (algorithm_name_ecmc, _, _, _, no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures_ecmc, _, _,
     no_of_jobs_ecmc, _, _) = run_script.get_config_data(base_config_file_ecmc)

    output_directory = sample_directory_4x4_metrop_low_temps.replace("/4x4_metrop_low_temps", "")
    sample_directories_metrop_low_temps = [f"{output_directory}/{length}x{length}_metrop_low_temps" for length in
                                           linear_system_sizes]
    sample_directories_metrop_lower_trans = [f"{output_directory}/{length}x{length}_metrop_lower_trans" for length in
                                             linear_system_sizes]
    sample_directories_metrop_upper_trans = [f"{output_directory}/{length}x{length}_metrop_upper_trans" for length in
                                             linear_system_sizes]
    sample_directories_metrop_high_temps = [f"{output_directory}/{length}x{length}_metrop_high_temps" for length in
                                            linear_system_sizes]
    sample_directories_ecmc = [f"{output_directory}/{length}x{length}_ecmc" for length in linear_system_sizes]
    pool = setup_pool(no_of_jobs_metrop_low_temps, max_no_of_cpus)

    temperatures_metrop = [*temperatures_metrop_low_temps, *temperatures_metrop_lower_trans,
                           *temperatures_metrop_upper_trans, *temperatures_metrop_high_temps]
    approx_transition_temperature = 0.887
    reduced_temperatures_metrop = [temperature / approx_transition_temperature for temperature in temperatures_metrop]
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
    inset_axis.set_xlim([0.89, 1.31]), inset_axis.set_ylim([2.0, 5.0 * 10 ** 3])
    inset_axis.set_yscale('log')
    [inset_axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]

    colors = ["black", "red", "blue", "green", "tab:brown", "magenta", "indigo"][:no_of_system_sizes]
    colors.reverse()

    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        print(f"Number of sites, N = {length}x{length}")
        """compute non-used physical time steps for our records"""
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_low_temps[system_size_index], temperatures_metrop_low_temps, length,
            no_of_observations_metrop_low_temps, no_of_jobs_metrop_low_temps, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_lower_trans[system_size_index], temperatures_metrop_lower_trans, length,
            no_of_observations_metrop_lower_trans, no_of_jobs_metrop_lower_trans, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_upper_trans[system_size_index], temperatures_metrop_upper_trans, length,
            no_of_observations_metrop_upper_trans, no_of_jobs_metrop_upper_trans, pool)
        compute_physical_time_steps(
            algorithm_name_metrop, external_global_moves_string, output_directory,
            sample_directories_metrop_high_temps[system_size_index], temperatures_metrop_high_temps, length,
            no_of_observations_metrop_high_temps, no_of_jobs_metrop_high_temps, pool)

        """compute non-used (in this script) twist probabilities for later use in make_twist_figs -- this optimises 
            our usage of the scratch space on BlueCrystal 4, ACRC, University of Bristol"""
        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_low_temps[system_size_index],
            temperatures_metrop_low_temps, length, no_of_observations_metrop_low_temps, no_of_jobs_metrop_low_temps,
            pool)
        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_lower_trans[system_size_index],
            temperatures_metrop_lower_trans, length, no_of_observations_metrop_lower_trans,
            no_of_jobs_metrop_lower_trans, pool)

        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_upper_trans[system_size_index],
            temperatures_metrop_upper_trans, length, no_of_observations_metrop_upper_trans,
            no_of_jobs_metrop_upper_trans, pool)

        _, _ = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_high_temps[system_size_index],
            temperatures_metrop_high_temps, length, no_of_observations_metrop_high_temps, no_of_jobs_metrop_high_temps,
            pool)

        cvms_metrop_low_temps, cvm_errors_metrop_low_temps = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_low_temps[system_size_index],
            no_of_equilibration_sweeps_metrop_low_temps, no_of_observations_metrop_low_temps,
            temperatures_metrop_low_temps, external_global_moves_string, no_of_jobs_metrop_low_temps, pool, length)
        cvms_metrop_lower_trans, cvm_errors_metrop_lower_trans = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_lower_trans[system_size_index],
            no_of_equilibration_sweeps_metrop_lower_trans, no_of_observations_metrop_lower_trans,
            temperatures_metrop_lower_trans, external_global_moves_string, no_of_jobs_metrop_lower_trans, pool, length)
        cvms_metrop_upper_trans, cvm_errors_metrop_upper_trans = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_upper_trans[system_size_index],
            no_of_equilibration_sweeps_metrop_upper_trans, no_of_observations_metrop_upper_trans,
            temperatures_metrop_upper_trans, external_global_moves_string, no_of_jobs_metrop_upper_trans, pool, length)
        cvms_metrop_high_temps, cvm_errors_metrop_high_temps = get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_high_temps[system_size_index],
            no_of_equilibration_sweeps_metrop_high_temps, no_of_observations_metrop_high_temps,
            temperatures_metrop_high_temps, external_global_moves_string, no_of_jobs_metrop_high_temps, pool, length)
        cvms_ecmc, cvm_errors_ecmc = get_cvm_mag_phase_vs_temperature(
            algorithm_name_ecmc, output_directory, sample_directories_ecmc[system_size_index],
            no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures_ecmc,
            external_global_moves_string, no_of_jobs_ecmc, pool, length)

        cvms_metrop = [*cvms_metrop_low_temps, *cvms_metrop_lower_trans, *cvms_metrop_upper_trans,
                       *cvms_metrop_high_temps]
        cvm_errors_metrop = [*cvm_errors_metrop_low_temps, *cvm_errors_metrop_lower_trans,
                             *cvm_errors_metrop_upper_trans, *cvm_errors_metrop_high_temps]

        axes[0].errorbar(reduced_temperatures_metrop, cvms_metrop, cvm_errors_metrop, marker=".", markersize=10,
                         color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length} Metrop.")
        axes[0].errorbar(reduced_temperatures_ecmc, cvms_ecmc, cvm_errors_ecmc, marker="*", markersize=8,
                         color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length} ECMC")
        inset_axis.errorbar(reduced_temperatures_metrop, cvms_metrop, cvm_errors_metrop, marker=".", markersize=10,
                            color=colors[system_size_index], linestyle="None")

    inverse_linear_system_sizes = [1.0 / length for length in linear_system_sizes]

    legend = axes[0].legend(loc="upper right", fontsize=9)
    legend.get_frame().set_edgecolor("k"), legend.get_frame().set_lw(3)
    fig.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_and_ecmc.pdf",
                bbox_inches="tight")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()


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
