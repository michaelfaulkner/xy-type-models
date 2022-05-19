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
cramer_von_mises = importlib.import_module("cramer_von_mises")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")
setup_scripts = importlib.import_module("setup_scripts")
sample_getter = importlib.import_module("sample_getter")
run_script = importlib.import_module("run")


def main():
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    # linear_system_sizes = [2 ** (index + 3) for index in range(5)]
    linear_system_sizes = [2 ** (index + 3) for index in range(1)]
    config_files_metrop = [f"config_files/cvm_figure/{value}x{value}_metrop.txt" for value in linear_system_sizes]
    config_files_ecmc = [f"config_files/cvm_figure/{value}x{value}_ecmc.txt" for value in linear_system_sizes]

    config_file_low_temp_all = f"config_files/cvm_figure/4x4_low_temp_all.txt"
    config_file_low_temp_local = f"config_files/cvm_figure/4x4_low_temp_local.txt"

    (algorithm_name_metrop, sample_directory_4x4_metrop, _, _, no_of_equilibration_sweeps_metrop,
     no_of_observations_metrop, temperatures, _, external_global_moves_string_all, no_of_jobs_metrop, _,
     max_no_of_cpus) = run_script.get_config_data(config_files_metrop[0])
    (algorithm_name_ecmc, _, _, _, no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, _, _, _, no_of_jobs_ecmc,
     _, _) = run_script.get_config_data(config_files_ecmc[0])

    (_, _, _, _, _, no_of_observations_metrop_low_temp, temperatures_low_temp, _, _, no_of_jobs_metrop_low_temp, _, _
     ) = run_script.get_config_data(config_file_low_temp_all)
    external_global_moves_string_local = run_script.get_config_data(config_file_low_temp_local)[8]

    output_directory = sample_directory_4x4_metrop.replace("/4x4_metrop", "")
    pool = setup_scripts.setup_pool(no_of_jobs_metrop, max_no_of_cpus)

    figure_cvm, axes_cvm = plt.subplots(1, 2, figsize=(10.0, 4.5))
    axes_cvm[0].text(1.1, 0.000000035, "(a)", fontsize=20)
    axes_cvm[1].text(0.3, -0.29, "(b)", fontsize=20)
    figure_cvm.tight_layout(w_pad=-0.5)
    axes_cvm[0].set_yscale('log')
    axes_cvm[0].set_xlabel(r"$(\beta J)^{- 1}$", fontsize=20, labelpad=3)
    axes_cvm[0].set_ylabel(r"$n \omega_{\phi_m,n}^2$", fontsize=20, labelpad=-4)
    axes_cvm[1].set_xlabel(r"$N^{- 1 / 2}$", fontsize=20, labelpad=3)
    axes_cvm[1].set_ylabel(r"$(\beta_{\rm int} \,\, J)^{-1}$", fontsize=20, labelpad=4)
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes_cvm]
    for axis_index, axis in enumerate(axes_cvm):
        axis.tick_params(which='both', direction='in', width=3)
        axis.tick_params(which='major', length=7, labelsize=18, pad=5)
        axis.tick_params(which='minor', length=4)

    additional_y_axis = axes_cvm[1].twinx()  # add a twinned y axis to the right-hand subplot
    additional_y_axis.set_ylabel(r"$\omega_{\phi_m,n}^{2,\rm{all}}(\beta J = 10) / "
                                 r"\omega_{\phi_m,n}^{2,\rm{local}}(\beta J = 10)$",
                                 fontsize=18, labelpad=4, color="red")
    additional_y_axis.tick_params(which='both', direction='in', width=3)
    additional_y_axis.tick_params(which='major', length=7, labelsize=18, pad=5)
    additional_y_axis.tick_params(which='minor', length=4)
    additional_y_axis.tick_params(axis='y', labelcolor='red')

    inset_axis = plt.axes([0.3, 0.66, 0.1525, 0.275])
    inset_axis.tick_params(which='both', direction='in', length=4, width=2, labelsize=12)
    inset_axis.set_xlim([0.86, 1.325])
    inset_axis.set_ylim([2.0, 5.0 * 10 ** 3])
    inset_axis.set_yscale('log')
    [inset_axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]

    figure_twist_probs, axis_twist_probs = plt.subplots(1, figsize=(5.0, 2.25))
    figure_twist_probs.tight_layout()
    axis_twist_probs.set_xlabel(r"$(\beta J)^{-1}$", fontsize=20, labelpad=3)
    axis_twist_probs.set_yscale('log')
    axis_twist_probs.set_ylim([9.0 * 10 ** (-8), 0.7])
    axis_twist_probs.set_yticks([10 ** (-7), 10 ** (-5), 10 ** (-3), 10 ** (-1)])
    axis_twist_probs.set_ylabel(r"$p(\rm{twist})$", fontsize=20, labelpad=4)
    [axis_twist_probs.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]
    axis_twist_probs.tick_params(which='both', direction='in', width=3)
    axis_twist_probs.tick_params(which='major', length=7, labelsize=18, pad=5)
    axis_twist_probs.tick_params(which='minor', length=4)

    colors = ["black", "red", "blue", "green", "yellow", "cyan"]

    cvm_ratios, cvm_ratio_errors = [], []
    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        try:
            with open(f"{output_directory}/probability_of_global_twists_{algorithm_name_metrop.replace('-', '_')}_"
                      f"{length}x{length}_sites_{no_of_observations_metrop}_obs_{no_of_jobs_metrop}_jobs.tsv",
                      "r") as output_file:
                output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                    if not line.startswith('#')]).transpose()
                twist_probabilities, twist_probability_errors = output_file_sans_header[1], output_file_sans_header[2]
        except IOError:
            twists_file = open(f"{output_directory}/probability_of_global_twists_"
                               f"{algorithm_name_metrop.replace('-', '_')}_{length}x{length}_sites_"
                               f"{no_of_observations_metrop}_obs_{no_of_jobs_metrop}_jobs.tsv", "w")
            twists_file.write("# temperature".ljust(30) + "p(twist)".ljust(30) + "p(twist) error" + "\n")
            twist_probabilities, twist_probability_errors = [], []
            for temperature_index, temperature in enumerate(temperatures):
                twist_probability_vs_job = np.array(pool.starmap(sample_getter.get_acceptance_rates, [
                    (f"{output_directory}/{length}x{length}_metrop/job_{job_index}", temperature_index)
                    for job_index in range(no_of_jobs_metrop)]))[:, 2]
                twist_probability, twist_probability_error = (np.mean(twist_probability_vs_job), np.std(
                    twist_probability_vs_job) / len(twist_probability_vs_job) ** 0.5)
                twist_probabilities.append(twist_probability)
                twist_probability_errors.append(twist_probability_error)
                twists_file.write(f"{temperature:.14e}".ljust(30) + f"{twist_probability:.14e}".ljust(30) +
                                  f"{twist_probability_error:.14e}" + "\n")
            twists_file.close()

        try:
            with open(f"{output_directory}/physical_time_steps_{algorithm_name_metrop.replace('-', '_')}_"
                      f"{external_global_moves_string_all}_{length}x{length}_sites_{no_of_observations_metrop}_obs_"
                      f"{no_of_jobs_metrop}_jobs.tsv", "r") as _:
                pass
        except IOError:
            physical_time_step_file = open(
                f"{output_directory}/physical_time_steps_{algorithm_name_metrop.replace('-', '_')}_"
                f"{external_global_moves_string_all}_{length}x{length}_sites_{no_of_observations_metrop}_obs_"
                f"{no_of_jobs_metrop}_jobs.tsv", "w")
            physical_time_step_file.write("# temperature".ljust(30) + "Delta t".ljust(30) + "Delta t error" + "\n")
            for temperature_index, temperature in enumerate(temperatures):
                physical_time_step_vs_job = pool.starmap(sample_getter.get_physical_time_step, [
                    (algorithm_name_metrop, f"{output_directory}/{length}x{length}_metrop/job_{job_index}",
                     temperature_index) for job_index in range(no_of_jobs_metrop)])
                physical_time_step_file.write(
                    f"{temperature:.14e}".ljust(30) + f"{np.mean(physical_time_step_vs_job):.14e}".ljust(30) +
                    f"{np.std(physical_time_step_vs_job) / len(physical_time_step_vs_job) ** 0.5:.14e}" + "\n")
            physical_time_step_file.close()

        cvm_metrops, cvm_metrop_errors = cramer_von_mises.get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, f"{output_directory}/{length}x{length}_metrop",
            no_of_equilibration_sweeps_metrop, no_of_observations_metrop, temperatures,
            external_global_moves_string_all, no_of_jobs_metrop, pool, length)
        cvm_ecmcs, cvm_ecmc_errors = cramer_von_mises.get_cvm_mag_phase_vs_temperature(
            algorithm_name_ecmc, output_directory, f"{output_directory}/{length}x{length}_ecmc",
            no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures,
            external_global_moves_string_all, no_of_jobs_ecmc, pool, length)
        low_temp_cvm_local, low_temp_cvm_local_error = cramer_von_mises.get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, f"{output_directory}/{length}x{length}_low_temp_local",
            no_of_equilibration_sweeps_metrop, no_of_observations_metrop, temperatures_low_temp,
            external_global_moves_string_local, no_of_jobs_metrop_low_temp, pool, length)
        low_temp_cvm_all, low_temp_cvm_all_error = cramer_von_mises.get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, f"{output_directory}/{length}x{length}_low_temp_all",
            no_of_equilibration_sweeps_metrop, no_of_observations_metrop, temperatures_low_temp,
            external_global_moves_string_all, no_of_jobs_metrop_low_temp, pool, length)

        cvm_ratios.append(low_temp_cvm_all[0] / low_temp_cvm_local[0])
        cvm_ratio_errors.append(math.sqrt(
            (low_temp_cvm_all_error[0] / low_temp_cvm_local[0]) ** 2 +
            (low_temp_cvm_all[0] * low_temp_cvm_local_error[0] / low_temp_cvm_local[0] ** 2) ** 2))

        axes_cvm[0].errorbar(temperatures, cvm_metrops, cvm_metrop_errors, marker=".", markersize=10,
                             color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length}")
        axes_cvm[0].errorbar(temperatures, cvm_ecmcs, cvm_ecmc_errors, marker="*", markersize=8,
                             color=colors[system_size_index], linestyle="None")
        inset_axis.errorbar(temperatures, cvm_metrops, cvm_metrop_errors, marker=".", markersize=8,
                            color=colors[system_size_index], linestyle="None")
        axis_twist_probs.errorbar(temperatures, twist_probabilities, twist_probability_errors, marker=".",
                                  markersize=10, color=colors[system_size_index], linestyle="None",
                                  label=fr"$N$ = {length}x{length}")

    '''cvm_ratio_file = open(
        f"{output_directory}/cvm_ratio_vs_system_size_{algorithm_name_metrop.replace('-', '_')}_temp_eq_"
        f"{temperatures_low_temp[0]:.4f}_{no_of_observations_metrop}_obs_{no_of_jobs_metrop}_jobs.tsv", "w")
    cvm_ratio_file.write("# no of sites".ljust(30) + "CvM ratio".ljust(30) + "CvM ratio error" + "\n")
    for linear_system_size_index, linear_system_size in enumerate(linear_system_sizes):
        cvm_ratio_file.write(f"{linear_system_size ** 2}".ljust(30) +
                             f"{cvm_ratios[linear_system_size_index]:.14e}".ljust(30) +
                             f"{cvm_ratio_errors[linear_system_size_index]:.14e}" + "\n")
    cvm_ratio_file.close()'''
    inverse_linear_system_sizes = [1.0 / length for length in linear_system_sizes]
    additional_y_axis.plot(inverse_linear_system_sizes, cvm_ratios, marker=".", markersize=8, color="red",
                           linestyle='None')

    legend_cvm = axes_cvm[0].legend(loc="center left", fontsize=10)
    legend_cvm.get_frame().set_edgecolor("k")
    legend_cvm.get_frame().set_lw(3)
    figure_cvm.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_and_ecmc"
                       f".pdf", bbox_inches="tight")

    legend_twist_probs = axis_twist_probs.legend(loc="lower right", fontsize=10)
    legend_twist_probs.get_frame().set_edgecolor("k")
    legend_twist_probs.get_frame().set_lw(3)
    figure_twist_probs.savefig(f"{output_directory}/probability_of_global_twists_xy_gaussian_noise_metropolis.pdf",
                               bbox_inches="tight")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: no positional arguments allowed.")
    else:
        main()
