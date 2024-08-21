from make_summary_statistic_vs_temperature_plot import get_means_and_errors
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
    base_config_file_electrolyte_all_moves = f"config_files/top_susc_figs/4x4_elec_all.txt"
    base_config_file_electrolyte_local_moves = f"config_files/top_susc_figs/4x4_elec_local.txt"
    base_config_file_hxy_all_moves = f"config_files/top_susc_figs/4x4_hxy_all.txt"
    base_config_file_xy_all_moves = f"config_files/top_susc_figs/4x4_xy_all.txt"
    (algorithm_name_electrolyte, sample_directory_4x4_electrolyte_all_moves, _, _, no_of_equilibration_sweeps,
     no_of_samples, temperatures_electrolyte, _, external_global_moves_string_all_moves, no_of_runs_electrolyte, _,
     max_no_of_cpus) = run_script.get_config_data(base_config_file_electrolyte_all_moves)
    (algorithm_name_hxy, _, _, _, _, _, temperatures_hxy, _, _, no_of_runs_hxy, _, _
     ) = run_script.get_config_data(base_config_file_hxy_all_moves)
    (algorithm_name_xy, _, _, _, _, _, temperatures_xy, _, _, no_of_runs_xy, _, _
     ) = run_script.get_config_data(base_config_file_xy_all_moves)
    external_global_moves_string_local_moves = run_script.get_config_data(base_config_file_electrolyte_local_moves)[8]

    algorithms = [algorithm_name_electrolyte, algorithm_name_hxy, algorithm_name_xy]
    no_of_runs = [no_of_runs_electrolyte, no_of_runs_hxy, no_of_runs_xy]
    """Define observables[i] as the necessary observables of algorithms[i]"""
    observables = [
        ["inverse_permittivity", "topological_susceptibility"],
        ["helicity_modulus", "global_defect_susceptibility", "hot_twist_relaxation_susceptibility"],
        ["helicity_modulus", "global_defect_susceptibility", "hot_twist_relaxation_susceptibility"]]
    observable_plotting_markers = [".", "*", "o"]
    system_size_plotting_colors = ["black", "red", "blue", "green", "magenta", "indigo"][:no_of_system_sizes]
    system_size_plotting_colors.reverse()

    output_directory = sample_directory_4x4_electrolyte_all_moves.replace("/4x4_elec_all", "")

    sample_directories_electrolyte_all_moves = [f"{output_directory}/{length}x{length}_elec_all" for length in
                                                linear_system_sizes]
    sample_directories_hxy_all_moves = [f"{output_directory}/{length}x{length}_hxy_all" for length in
                                        linear_system_sizes]
    sample_directories_xy_all_moves = [f"{output_directory}/{length}x{length}_xy_all" for length in linear_system_sizes]
    sample_directories_by_algo_all_moves = [sample_directories_electrolyte_all_moves, sample_directories_hxy_all_moves,
                                            sample_directories_xy_all_moves]

    sample_directories_electrolyte_local_moves = [f"{output_directory}/{length}x{length}_elec_local" for length in
                                                  linear_system_sizes]
    sample_directories_hxy_local_moves = [f"{output_directory}/{length}x{length}_hxy_local" for length in
                                          linear_system_sizes]
    sample_directories_xy_local_moves = [f"{output_directory}/{length}x{length}_xy_local" for length in
                                         linear_system_sizes]
    sample_directories_by_algo_local_moves = [sample_directories_electrolyte_local_moves,
                                              sample_directories_hxy_local_moves, sample_directories_xy_local_moves]

    model_temperatures = [temperatures_electrolyte, temperatures_hxy, temperatures_xy]
    approx_transition_temperatures = [1.351, 1.351, 0.887]
    reduced_model_temperatures = [
        [temperature / approx_transition_temperatures[0] for temperature in temperatures_electrolyte],
        [temperature / approx_transition_temperatures[1] for temperature in temperatures_hxy],
        [temperature / approx_transition_temperatures[2] for temperature in temperatures_xy]]

    start_time = time.time()
    make_topological_susceptibility_plots(algorithms, observables, model_temperatures, reduced_model_temperatures,
                                          linear_system_sizes, no_of_samples, no_of_runs, output_directory,
                                          sample_directories_by_algo_all_moves, external_global_moves_string_all_moves,
                                          observable_plotting_markers, system_size_plotting_colors)
    make_topological_stability_plots(algorithms, observables, model_temperatures, reduced_model_temperatures,
                                     linear_system_sizes, no_of_samples, no_of_runs, output_directory,
                                     sample_directories_by_algo_local_moves, sample_directories_by_algo_all_moves,
                                     external_global_moves_string_local_moves, external_global_moves_string_all_moves,
                                     observable_plotting_markers, system_size_plotting_colors)
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


def make_topological_susceptibility_plots(algorithms, observables, model_temperatures, reduced_model_temperatures,
                                          linear_system_sizes, no_of_samples, no_of_runs, output_directory,
                                          sample_directories_by_algo, external_global_moves_string_all_moves,
                                          observable_plotting_markers, system_size_plotting_colors):
    fig, axes = make_three_empty_subfigures(r"$\beta L^2 {\rm Var}[{\rm X}] / (2J)$")
    for algorithm_index, algorithm in enumerate(algorithms):
        for system_size_index, length in enumerate(linear_system_sizes):
            no_of_sites_string = f"{length}x{length}_sites"
            """first compute mean acceptance rates across the runs"""
            get_means_and_errors("acceptance_rates", algorithm, model_temperatures[algorithm_index],
                                 no_of_sites_string, no_of_samples, external_global_moves_string_all_moves,
                                 output_directory, sample_directories_by_algo[algorithm_index][system_size_index],
                                 no_of_runs[algorithm_index])
            for observable_index, observable in enumerate(observables[algorithm_index]):
                means, errors = get_means_and_errors(observable, algorithm, model_temperatures[algorithm_index],
                                                     no_of_sites_string, no_of_samples,
                                                     external_global_moves_string_all_moves, output_directory,
                                                     sample_directories_by_algo[algorithm_index][system_size_index],
                                                     no_of_runs[algorithm_index])
                if observable_index > 0:
                    axes[algorithm_index].errorbar(reduced_model_temperatures[algorithm_index], means, errors,
                                                   marker=observable_plotting_markers[observable_index - 1],
                                                   markersize=10, color=system_size_plotting_colors[system_size_index],
                                                   linestyle="dashed")
    fig.savefig(f"{output_directory}/topological_susceptibilities_with_global_moves.pdf", bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


def make_topological_stability_plots(algorithms, observables, model_temperatures, reduced_model_temperatures,
                                     linear_system_sizes, no_of_samples, no_of_runs, output_directory,
                                     sample_directories_by_algo_local_moves, sample_directories_by_algo_all_moves,
                                     external_global_moves_string_local_moves, external_global_moves_string_all_moves,
                                     observable_plotting_markers, system_size_plotting_colors):
    fig, axes = make_three_empty_subfigures(r"$1 - \sqrt{\langle s_{\rm X}^2 \rangle / {\rm Var}[{\rm X}]}$")
    for algorithm_index, algorithm in enumerate(algorithms):
        for system_size_index, length in enumerate(linear_system_sizes):
            no_of_sites_string = f"{length}x{length}_sites"
            """first compute mean acceptance rates across the runs"""
            get_means_and_errors("acceptance_rates", algorithm, model_temperatures[algorithm_index], no_of_sites_string,
                                 no_of_samples, external_global_moves_string_local_moves, output_directory,
                                 sample_directories_by_algo_local_moves[algorithm_index][system_size_index],
                                 no_of_runs[algorithm_index])
            get_means_and_errors("acceptance_rates", algorithm, model_temperatures[algorithm_index], no_of_sites_string,
                                 no_of_samples, external_global_moves_string_all_moves, output_directory,
                                 sample_directories_by_algo_all_moves[algorithm_index][system_size_index],
                                 no_of_runs[algorithm_index])

            for observable_index, observable in enumerate(observables[algorithm_index]):
                output_file_string = (f"{output_directory}/{observable}_ratio_{algorithm.replace('-', '_')}_"
                                      f"{no_of_sites_string}_{no_of_runs[algorithm_index]}x{no_of_samples}_obs.tsv")
                try:
                    with open(output_file_string, "r") as output_file:
                        output_file_sans_header = np.array([
                            np.fromstring(line, dtype=float, sep='\t')
                            for line in output_file if not line.startswith('#')]).transpose()
                        chi_ratios = output_file_sans_header[1]
                except IOError:
                    sample_means_local_moves, sample_errors_local_moves = get_means_and_errors(
                        observable, algorithm, model_temperatures[algorithm_index], no_of_sites_string, no_of_samples,
                        external_global_moves_string_local_moves, output_directory,
                        sample_directories_by_algo_local_moves[algorithm_index][system_size_index],
                        no_of_runs[algorithm_index])
                    sample_means_all_moves, sample_errors_all_moves = get_means_and_errors(
                        observable, algorithm, model_temperatures[algorithm_index], no_of_sites_string, no_of_samples,
                        external_global_moves_string_all_moves, output_directory,
                        sample_directories_by_algo_all_moves[algorithm_index][system_size_index],
                        no_of_runs[algorithm_index])
                    chi_ratios = [sample_mean_local_moves / sample_means_all_moves[temperature_index]
                                  for temperature_index, sample_mean_local_moves in enumerate(sample_means_local_moves)]
                    chi_ratio_errors = [math.sqrt(
                        (sample_errors_local_moves[temperature_index] / sample_means_all_moves[temperature_index]) ** 2
                        + (sample_mean_local_moves * sample_errors_all_moves[temperature_index] /
                           sample_means_all_moves[temperature_index] ** 2) ** 2)
                        for temperature_index, sample_mean_local_moves in enumerate(sample_means_local_moves)]

                    output_file = open(output_file_string, "w")
                    output_file.write("# *** NB, chi represents the observable in the file name ***" + "\n")
                    output_file.write("# temperature".ljust(30) + "chi ratio".ljust(30) + "chi ratio error".ljust(30) +
                                      "chi (local only)".ljust(30) + "chi error (local only)".ljust(30) +
                                      "chi (all moves)".ljust(30) + "chi error (all moves)".ljust(30) + "\n")
                    for temperature_index, temperature in enumerate(model_temperatures[algorithm_index]):
                        output_file.write(f"{temperature:.14e}".ljust(30) +
                                          f"{chi_ratios[temperature_index]:.14e}".ljust(30) +
                                          f"{chi_ratio_errors[temperature_index]:.14e}".ljust(30) +
                                          f"{sample_means_local_moves[temperature_index]:.14e}".ljust(30) +
                                          f"{sample_errors_local_moves[temperature_index]:.14e}".ljust(30) +
                                          f"{sample_means_all_moves[temperature_index]:.14e}".ljust(30) +
                                          f"{sample_errors_all_moves[temperature_index]:.14e}".ljust(30) + "\n")
                    output_file.close()

                topological_stabilities = [1.0 - chi_ratio ** 0.5 for chi_ratio in chi_ratios]
                axes[algorithm_index].plot(reduced_model_temperatures[algorithm_index], topological_stabilities,
                                           marker=observable_plotting_markers[observable_index], markersize=10,
                                           color=system_size_plotting_colors[system_size_index], linestyle="dashed")

    fig.savefig(f"{output_directory}/topological_stabilities.pdf", bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


def make_three_empty_subfigures(y_axis_label):
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))
    fig.tight_layout(w_pad=2.5)
    axes[0].set_ylabel(y_axis_label, fontsize=25, labelpad=-1.0)
    for axis_index, axis in enumerate(axes):
        axis.tick_params(which='both', direction='in', width=3)
        axis.tick_params(which='major', length=7, labelsize=22.5, pad=3)
        axis.tick_params(which='minor', length=4)
        [axis.spines[spine].set_linewidth(3.75) for spine in ["top", "bottom", "left", "right"]]
        if axis_index < 2:
            axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=20, labelpad=-0.5)
        else:
            axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta$", fontsize=20, labelpad=-0.5)
        """
        axis.set_ylim(0.0, 1.0)
        axis.yaxis.set_major_locator(ticker.MultipleLocator(base=0.5))
        axis.yaxis.set_major_formatter('{x:.1e}')
        axis.axes.yaxis.set_ticklabels([])
        """
    return fig, axes


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("InterfaceError: At most one positional arguments permitted.  None are required but you may "
                        "provide no_of_system_sizes, which must be an integer greater than 0 and less than 7 (default "
                        "value is 6).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be no_of_system_sizes - which must be an integer greater "
              "than 0 and less than 7 (default value is 6).")
        chosen_no_of_system_sizes = int(sys.argv[1])
        if chosen_no_of_system_sizes < 1 or chosen_no_of_system_sizes > 5:
            raise Exception("InterfaceError: no_of_system_sizes must be an integer greater than 0 and less than 7 "
                            "(default value is 6).")
        main(chosen_no_of_system_sizes)
    else:
        print("No positional arguments provided.  None are required but you may provide no_of_system_sizes, which must "
              "be an integer greater than 0 and less than 7 (default value is 6).")
        main()
