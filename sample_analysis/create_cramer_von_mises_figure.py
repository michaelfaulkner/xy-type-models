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


def main(no_of_system_sizes=7):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    linear_system_sizes = [2 ** (index + 2) for index in range(no_of_system_sizes)]
    base_config_file_metrop_high_temps = f"config_files/cvm_figure/4x4_metrop.txt"
    base_config_file_ecmc = f"config_files/cvm_figure/4x4_ecmc.txt"
    base_config_file_low_temp_all = f"config_files/cvm_figure/4x4_low_temp_all.txt"
    base_config_file_low_temp_local = f"config_files/cvm_figure/4x4_low_temp_local.txt"

    """n.b., temperatures = [*temperatures_low_temp, *temperatures_high_temps] in the following lines"""
    (algorithm_name_metrop, sample_directory_4x4_metrop_high_temps, _, _, no_of_equilibration_sweeps_metrop,
     no_of_observations_metrop, temperatures_high_temps, _, external_global_moves_string_all,
     no_of_jobs_metrop_high_temps, _, max_no_of_cpus) = run_script.get_config_data(base_config_file_metrop_high_temps)
    (algorithm_name_ecmc, _, _, _, no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures, _, _,
     no_of_jobs_ecmc, _, _) = run_script.get_config_data(base_config_file_ecmc)
    (_, _, _, _, _, _, temperatures_low_temp, _, _, no_of_jobs_metrop_low_temp, _, _
     ) = run_script.get_config_data(base_config_file_low_temp_all)
    external_global_moves_string_local = run_script.get_config_data(base_config_file_low_temp_local)[8]

    output_directory = sample_directory_4x4_metrop_high_temps.replace("/4x4_metrop", "")
    sample_directories_metrop_high_temps = [f"{output_directory}/{length}x{length}_metrop" for length in
                                            linear_system_sizes]
    sample_directories_low_temp_local = [f"{output_directory}/{length}x{length}_low_temp_local" for length in
                                         linear_system_sizes]
    sample_directories_low_temp_all = [f"{output_directory}/{length}x{length}_low_temp_all" for length in
                                       linear_system_sizes]
    sample_directories_ecmc = [f"{output_directory}/{length}x{length}_ecmc" for length in linear_system_sizes]
    pool = setup_scripts.setup_pool(no_of_jobs_metrop_low_temp, max_no_of_cpus)

    figure_cvm, axes_cvm = plt.subplots(1, 2, figsize=(10.0, 4.5), gridspec_kw={'width_ratios': [1.2, 1.0]})
    figure_cvm.text(0.0025, 0.925, "d", fontsize=20, weight='bold')
    figure_cvm.text(0.52, 0.925, "e", fontsize=20, weight='bold')
    figure_cvm.tight_layout(w_pad=3.0)
    axes_cvm[0].set_yscale('log')
    axes_cvm[0].set_xlabel(r"$1 / (\beta J)$", fontsize=20, labelpad=3)
    axes_cvm[0].set_ylabel(r"$n \omega_{\phi_m,n}^2$", fontsize=20, labelpad=-4)
    axes_cvm[1].set_xlabel(r"$1 / N$", fontsize=20, labelpad=3)
    axes_cvm[1].set_ylabel(r"$1 / (\beta_{\rm int} \,\, J)$", fontsize=20, labelpad=1)
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes_cvm]
    for axis_index, axis in enumerate(axes_cvm):
        axis.tick_params(which='both', direction='in', width=3)
        axis.tick_params(which='major', length=7, labelsize=18, pad=5)
        axis.tick_params(which='minor', length=4)

    additional_y_axis = axes_cvm[1].twinx()  # add a twinned y axis to the right-hand subplot
    additional_y_axis.set_ylabel(r"$\chi_(\beta J = 10)$", fontsize=20, labelpad=1, color="red")
    additional_y_axis.tick_params(which='both', direction='in', width=3)
    additional_y_axis.tick_params(which='major', length=7, labelsize=18, pad=5)
    additional_y_axis.tick_params(which='minor', length=4)
    additional_y_axis.tick_params(axis='y', labelcolor='red')

    inset_axis = plt.axes([0.345, 0.66, 0.1525, 0.275])
    inset_axis.tick_params(which='both', direction='in', length=4, width=2, labelsize=12)
    inset_axis.set_xlim([0.86, 1.325])
    inset_axis.set_ylim([2.0, 5.0 * 10 ** 3])
    inset_axis.set_yscale('log')
    [inset_axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]

    figure_twist_probs, axis_twist_probs = plt.subplots(1, figsize=(5.0, 2.25))
    figure_twist_probs.tight_layout()
    axis_twist_probs.set_xlabel(r"$1 / (\beta J)$", fontsize=20, labelpad=3)
    axis_twist_probs.set_yscale('log')
    axis_twist_probs.set_ylim([9.0 * 10 ** (-8), 0.7])
    axis_twist_probs.set_yticks([10 ** (-7), 10 ** (-5), 10 ** (-3), 10 ** (-1)])
    axis_twist_probs.set_ylabel(r"$p(\rm{twist})$", fontsize=20, labelpad=4)
    [axis_twist_probs.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]
    axis_twist_probs.tick_params(which='both', direction='in', width=3)
    axis_twist_probs.tick_params(which='major', length=7, labelsize=18, pad=5)
    axis_twist_probs.tick_params(which='minor', length=4)

    colors = ["magenta", "cyan", "yellow", "green", "blue", "red", "black"]

    cvm_ratios, cvm_ratio_errors = [], []
    sample_variance_vs_system_size_low_temp_local = []
    sample_variance_error_vs_system_size_low_temp_local = []
    sample_variance_vs_system_size_low_temp_all = []
    sample_variance_error_vs_system_size_low_temp_all = []
    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        print(f"No of sites, N = {length}x{length}")

        twist_probabilities_low_temp, twist_probability_errors_low_temp = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_low_temp_all[system_size_index],
            temperatures_low_temp, length, no_of_observations_metrop, no_of_jobs_metrop_low_temp, pool)
        twist_probabilities_high_temps, twist_probability_errors_high_temps = get_twist_probabilities_and_errors(
            algorithm_name_metrop, output_directory, sample_directories_metrop_high_temps[system_size_index],
            temperatures_high_temps, length, no_of_observations_metrop, no_of_jobs_metrop_high_temps, pool)
        twist_probabilities = [*twist_probabilities_low_temp, *twist_probabilities_high_temps]
        twist_probability_errors = [*twist_probability_errors_low_temp, *twist_probability_errors_high_temps]

        compute_physical_time_steps(algorithm_name_metrop, external_global_moves_string_all, output_directory,
                                    sample_directories_low_temp_all[system_size_index], temperatures_low_temp, length,
                                    no_of_observations_metrop, no_of_jobs_metrop_low_temp, pool)
        compute_physical_time_steps(algorithm_name_metrop, external_global_moves_string_all, output_directory,
                                    sample_directories_metrop_high_temps[system_size_index], temperatures_high_temps,
                                    length, no_of_observations_metrop, no_of_jobs_metrop_high_temps, pool)

        cvms_metrop_high_temps, cvm_errors_metrop_high_temps = cramer_von_mises.get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_metrop_high_temps[system_size_index],
            no_of_equilibration_sweeps_metrop, no_of_observations_metrop, temperatures_high_temps,
            external_global_moves_string_all, no_of_jobs_metrop_high_temps, pool, length)
        cvms_ecmc, cvm_errors_ecmc = cramer_von_mises.get_cvm_mag_phase_vs_temperature(
            algorithm_name_ecmc, output_directory, sample_directories_ecmc[system_size_index],
            no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures,
            external_global_moves_string_all, no_of_jobs_ecmc, pool, length)
        low_temp_cvm_local, low_temp_cvm_local_error = cramer_von_mises.get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_low_temp_local[system_size_index],
            no_of_equilibration_sweeps_metrop, no_of_observations_metrop, temperatures_low_temp,
            external_global_moves_string_local, no_of_jobs_metrop_low_temp, pool, length)
        low_temp_cvm_all, low_temp_cvm_all_error = cramer_von_mises.get_cvm_mag_phase_vs_temperature(
            algorithm_name_metrop, output_directory, sample_directories_low_temp_all[system_size_index],
            no_of_equilibration_sweeps_metrop, no_of_observations_metrop, temperatures_low_temp,
            external_global_moves_string_all, no_of_jobs_metrop_low_temp, pool, length)

        cvms_metrop, cvm_errors_metrop = [*low_temp_cvm_all, *cvms_metrop_high_temps], [*low_temp_cvm_all_error,
                                                                                        *cvm_errors_metrop_high_temps]
        cvm_ratios.append(low_temp_cvm_all[0] / low_temp_cvm_local[0])
        cvm_ratio_errors.append(math.sqrt(
            (low_temp_cvm_all_error[0] / low_temp_cvm_local[0]) ** 2 +
            (low_temp_cvm_all[0] * low_temp_cvm_local_error[0] / low_temp_cvm_local[0] ** 2) ** 2))

        (sample_variances_low_temp_local,
         sample_variance_errors_low_temp_local) = get_mag_phase_sample_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string_local, output_directory,
            sample_directories_low_temp_local[system_size_index], temperatures_low_temp, length,
            no_of_equilibration_sweeps_metrop, no_of_observations_metrop, no_of_jobs_metrop_low_temp, pool)
        sample_variances_low_temp_all, sample_variance_errors_low_temp_all = get_mag_phase_sample_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string_all, output_directory,
            sample_directories_low_temp_all[system_size_index], temperatures_low_temp, length,
            no_of_equilibration_sweeps_metrop, no_of_observations_metrop, no_of_jobs_metrop_low_temp, pool)
        _, _ = get_mag_phase_sample_variances_and_errors(
            algorithm_name_metrop, external_global_moves_string_all, output_directory,
            sample_directories_metrop_high_temps[system_size_index], temperatures_high_temps, length,
            no_of_equilibration_sweeps_metrop, no_of_observations_metrop, no_of_jobs_metrop_high_temps, pool)
        _, _ = get_mag_phase_sample_variances_and_errors(
            algorithm_name_ecmc, external_global_moves_string_all, output_directory,
            sample_directories_ecmc[system_size_index], temperatures, length, no_of_equilibration_sweeps_ecmc,
            no_of_observations_ecmc, no_of_jobs_ecmc, pool)

        sample_variance_vs_system_size_low_temp_local.append(sample_variances_low_temp_local[0])
        sample_variance_error_vs_system_size_low_temp_local.append(sample_variance_errors_low_temp_local[0])
        sample_variance_vs_system_size_low_temp_all.append(sample_variances_low_temp_all[0])
        sample_variance_error_vs_system_size_low_temp_all.append(sample_variance_errors_low_temp_all[0])

        axes_cvm[0].errorbar(temperatures, cvms_metrop, cvm_errors_metrop, marker=".", markersize=10,
                             color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length}")
        axes_cvm[0].errorbar(temperatures, cvms_ecmc, cvm_errors_ecmc, marker="*", markersize=8,
                             color=colors[system_size_index], linestyle="None")
        inset_axis.errorbar(temperatures, cvms_metrop, cvm_errors_metrop, marker=".", markersize=8,
                            color=colors[system_size_index], linestyle="None")
        axis_twist_probs.errorbar(temperatures, twist_probabilities, twist_probability_errors, marker=".",
                                  markersize=10, color=colors[system_size_index], linestyle="None",
                                  label=fr"$N$ = {length}x{length}")

    cvm_ratio_file = open(
        f"{output_directory}/cvm_ratio_vs_system_size_{algorithm_name_metrop.replace('-', '_')}_temp_eq_"
        f"{temperatures_low_temp[0]:.4f}_{no_of_observations_metrop}_obs_{no_of_jobs_metrop_low_temp}_jobs.tsv", "w")
    cvm_ratio_file.write("# no of sites".ljust(30) + "CvM ratio".ljust(30) + "CvM ratio error" + "\n")
    for linear_system_size_index, linear_system_size in enumerate(linear_system_sizes):
        cvm_ratio_file.write(f"{linear_system_size ** 2}".ljust(30) +
                             f"{cvm_ratios[linear_system_size_index]:.14e}".ljust(30) +
                             f"{cvm_ratio_errors[linear_system_size_index]:.14e}" + "\n")
    cvm_ratio_file.close()

    inverse_system_sizes = [1.0 / length ** 2 for length in linear_system_sizes]

    # TODO add errors in the following plotting functions
    additional_y_axis.plot(inverse_system_sizes, sample_variance_vs_system_size_low_temp_local, marker=".",
                           markersize=8, color="red", linestyle='None',
                           label=r"$\chi(\beta J) = s_{\phi_m,n = 10^6}^2(\beta J)$; local")
    additional_y_axis.plot(inverse_system_sizes, sample_variance_vs_system_size_low_temp_all, marker="d", markersize=8,
                           color="red", linestyle='None',
                           label=r"$\chi(\beta J) = s_{\phi_m,n = 10^6}^2(\beta J)$; all")
    additional_y_axis.plot(inverse_system_sizes, cvm_ratios, marker="*", markersize=8, color="red",
                           linestyle='None', label=r"$\chi(\beta J) = \omega_{\phi_m,n}^{2,\rm{all}}(\beta J) / "
                                                   r"\omega_{\phi_m,n}^{2,\rm{local}}(\beta J)$")

    legends_cvm = [axes_cvm[0].legend(loc="center left", fontsize=10),
                   additional_y_axis.legend(loc="upper left", fontsize=10)]
    [legend.get_frame().set_edgecolor("k") for legend in legends_cvm]
    [legend.get_frame().set_lw(3) for legend in legends_cvm]
    figure_cvm.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_and_ecmc"
                       f".pdf", bbox_inches="tight")

    legend_twist_probs = axis_twist_probs.legend(loc="lower right", fontsize=10)
    legend_twist_probs.get_frame().set_edgecolor("k")
    legend_twist_probs.get_frame().set_lw(3)
    figure_twist_probs.savefig(f"{output_directory}/probability_of_global_twists_xy_gaussian_noise_metropolis.pdf",
                               bbox_inches="tight")

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
            twist_probability_vs_job = np.array(pool.starmap(sample_getter.get_acceptance_rates, [
                (f"{sample_directory}/job_{job_index}", temperature_index) for job_index in range(no_of_jobs)]))[:, 2]
            twist_probability, twist_probability_error = (np.mean(twist_probability_vs_job), np.std(
                twist_probability_vs_job) / len(twist_probability_vs_job) ** 0.5)
            twist_probabilities.append(twist_probability)
            twist_probability_errors.append(twist_probability_error)
            twists_file.write(f"{temperature:.14e}".ljust(30) + f"{twist_probability:.14e}".ljust(30) +
                              f"{twist_probability_error:.14e}" + "\n")
        twists_file.close()
    return twist_probabilities, twist_probability_errors


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
            physical_time_step_vs_job = pool.starmap(sample_getter.get_physical_time_step, [
                (algorithm_name, f"{sample_directory}/job_{job_index}", temperature_index) for job_index in
                range(no_of_jobs)])
            physical_time_step_file.write(
                f"{temperature:.14e}".ljust(30) + f"{np.mean(physical_time_step_vs_job):.14e}".ljust(30) +
                f"{np.std(physical_time_step_vs_job) / len(physical_time_step_vs_job) ** 0.5:.14e}" + "\n")
        physical_time_step_file.close()


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
    return np.var(sample_getter.get_magnetisation_phase(sample_directory, temperature, temperature_index,
                                                        no_of_sites)[no_of_equilibration_sweeps:])


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("InterfaceError: At most one positional arguments permitted.  None are required but you may "
                        "provide no_of_system_sizes, which must be an integer greater than 0 and less than 8 (default "
                        "value is 7).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be no_of_system_sizes - which must be an integer greater "
              "than 0 and less than 8 (default value is 7).")
        chosen_no_of_system_sizes = int(sys.argv[1])
        if chosen_no_of_system_sizes < 1 or chosen_no_of_system_sizes > 7:
            raise Exception("InterfaceError: no_of_system_sizes must be an integer greater than 0 and less than 8 "
                            "(default value is 7).")
        main(chosen_no_of_system_sizes)
    else:
        print("No positional arguments orivuded.  None are required but you may provide no_of_system_sizes, which must "
              "be an integer greater than 0 and less than 8 (default value is 7).")
        main()
