import importlib
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
    linear_system_sizes = [2 ** (index + 3) for index in range(4)]
    config_files_metrop = [f"config_files/cvm_figure/{value}x{value}_metrop.txt" for value in linear_system_sizes]
    config_files_ecmc = [f"config_files/cvm_figure/{value}x{value}_ecmc.txt" for value in linear_system_sizes]

    (algorithm_name_metrop, sample_directory_8x8_metrop, _, _, no_of_equilibration_sweeps_metrop,
     no_of_observations_metrop, temperatures, _, external_global_moves_string, no_of_jobs_metrop, _,
     max_no_of_cpus) = run_script.get_config_data(config_files_metrop[0])
    (algorithm_name_ecmc, _, _, _, no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, _, _, _, no_of_jobs_ecmc,
     _, _) = run_script.get_config_data(config_files_ecmc[0])
    output_directory = sample_directory_8x8_metrop.replace("/8x8_metrop", "")
    pool = setup_scripts.setup_pool(no_of_jobs_metrop, max_no_of_cpus)

    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.5))
    axes[0].text(1.1, 0.000000035, "(a)", fontsize=20)
    axes[1].text(0.45, -0.29, "(b)", fontsize=20)
    figure.tight_layout(w_pad=-0.5)
    axes[0].set_yscale('log')
    axes[0].set_xlabel(r"$1 / (\beta J)$", fontsize=20, labelpad=3)
    axes[0].set_ylabel(r"$n \omega_{\phi_m,n}^2$", fontsize=20, labelpad=-4)
    axes[1].set_xlabel(r"$1 / N$", fontsize=20, labelpad=3)
    axes[1].set_ylabel(r"$1 / (\beta_{\rm int} \,\, J)$", fontsize=20, labelpad=4)
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes]
    for axis_index, axis in enumerate(axes):
        axis.tick_params(which='both', direction='in', width=3)
        axis.tick_params(which='major', length=7, labelsize=18, pad=5)
        axis.tick_params(which='minor', length=4)

    inset_axis_1 = plt.axes([0.3, 0.66, 0.1525, 0.275])
    inset_axis_1.tick_params(which='both', direction='in', length=4, width=2, labelsize=12)
    inset_axis_1.set_xlim([0.86, 1.325])
    inset_axis_1.set_ylim([2.0, 5.0 * 10 ** 3])
    inset_axis_1.set_yscale('log')
    [inset_axis_1.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]

    inset_axis_2 = plt.axes([0.8, 0.2925, 0.165, 0.26])
    inset_axis_2.tick_params(which='both', direction='in', length=4, width=2, labelsize=12)
    inset_axis_2.xaxis.set_label_position("top")
    inset_axis_2.xaxis.tick_top()
    inset_axis_2.set_xticks([0.0, 1.0, 2.0])
    inset_axis_2.set_xlabel(r"$1 / (\beta J)$", fontsize=12, labelpad=2)
    inset_axis_2.set_yscale('log')
    inset_axis_2.set_yticks([10 ** (-5), 10 ** (-3), 10 ** (-1)])
    inset_axis_2.set_ylabel(r"$p(\rm{twist})$", fontsize=12, labelpad=2)
    [inset_axis_2.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]

    colors = ["black", "red", "blue", "green", "yellow", "cyan"]

    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        try:
            with open(f"{output_directory}/probability_of_global_twists_{algorithm_name_metrop.replace('-', '_')}_"
                      f"{length}x{length}_sites.tsv", "r") as output_file:
                output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                    if not line.startswith('#')]).transpose()
                twist_probabilities, twist_probability_errors = output_file_sans_header[1], output_file_sans_header[2]
        except IOError:
            twists_file = open(f"{output_directory}/probability_of_global_twists_"
                               f"{algorithm_name_metrop.replace('-', '_')}_{length}x{length}_sites.tsv", "w")
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
                      f"{external_global_moves_string}_{length}x{length}_sites.tsv", "r") as _:
                pass
        except IOError:
            physical_time_step_file = open(
                f"{output_directory}/physical_time_steps_{algorithm_name_metrop.replace('-', '_')}_"
                f"{external_global_moves_string}_{length}x{length}_sites.tsv", "w")
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
            external_global_moves_string, no_of_jobs_metrop, pool, length)
        cvm_ecmcs, cvm_ecmc_errors = cramer_von_mises.get_cvm_mag_phase_vs_temperature(
            algorithm_name_ecmc, output_directory, f"{output_directory}/{length}x{length}_ecmc",
            no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc, temperatures,
            external_global_moves_string, no_of_jobs_ecmc, pool, length)
        axes[0].errorbar(temperatures, cvm_metrops, cvm_metrop_errors, marker=".", markersize=10,
                         color=colors[system_size_index], linestyle="None", label=fr"$N$ = {length}x{length}")
        axes[0].errorbar(temperatures, cvm_ecmcs, cvm_ecmc_errors, marker="*", markersize=8,
                         color=colors[system_size_index], linestyle="None")
        inset_axis_1.errorbar(temperatures, cvm_metrops, cvm_metrop_errors, marker=".", markersize=8,
                              color=colors[system_size_index], linestyle="None")
        inset_axis_2.errorbar(temperatures, twist_probabilities, twist_probability_errors, marker=".",
                              markersize=8, color=colors[system_size_index], linestyle="None")

    legend = axes[0].legend(loc="center left", fontsize=10)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(3)
    figure.savefig(f"{output_directory}/magnetisation_phase_cramervonmises_xy_gaussian_noise_metropolis_and_ecmc.pdf",
                   bbox_inches="tight")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: no positional arguments allowed.")
    else:
        main()
