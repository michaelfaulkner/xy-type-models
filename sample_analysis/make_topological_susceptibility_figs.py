from make_summary_statistic_vs_temperature_plot import get_means_and_errors
from matplotlib.legend_handler import HandlerLine2D
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
    observables = [["zero_mode_susceptibility", "topological_susceptibility"],
                   ["josephson_susceptibility", "global_defect_susceptibility", "hot_twist_relaxation_susceptibility"],
                   ["josephson_susceptibility", "global_defect_susceptibility", "hot_twist_relaxation_susceptibility",
                    "emergent_field_zero_mode_susceptibility"]]
    observable_plotting_markers = [".", "*"]
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
    make_helicity_plots(algorithms, model_temperatures, reduced_model_temperatures, linear_system_sizes,
                        no_of_samples, no_of_runs, output_directory, sample_directories_by_algo_all_moves,
                        external_global_moves_string_all_moves, system_size_plotting_colors)
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


def make_helicity_plots(algorithms, model_temperatures, reduced_model_temperatures, linear_system_sizes,
                        no_of_samples, no_of_runs, output_directory, sample_directories_by_algo,
                        external_global_moves_string_all_moves, system_size_plotting_colors):
    fig, axes = make_three_empty_subfigures(r"$\Gamma / \, J$", 25.0)
    [axes[axis_index].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) for axis_index in [0, 1]]
    axes[2].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    observables = ["inverse_permittivity", "helicity_modulus", "helicity_modulus"]
    for algorithm_index, algorithm in enumerate(algorithms):
        for system_size_index, length in enumerate(linear_system_sizes):
            no_of_sites_string = f"{length}x{length}_sites"
            """first compute mean acceptance rates across the runs"""
            get_means_and_errors("acceptance_rates", algorithm, model_temperatures[algorithm_index],
                                 no_of_sites_string, no_of_samples, external_global_moves_string_all_moves,
                                 output_directory, sample_directories_by_algo[algorithm_index][system_size_index],
                                 no_of_runs[algorithm_index])
            means, errors = get_means_and_errors(
                observables[algorithm_index], algorithm, model_temperatures[algorithm_index], no_of_sites_string,
                no_of_samples, external_global_moves_string_all_moves, output_directory,
                sample_directories_by_algo[algorithm_index][system_size_index], no_of_runs[algorithm_index])
            axes[algorithm_index].errorbar(reduced_model_temperatures[algorithm_index], means, errors, marker=".",
                                           markersize=12, color=system_size_plotting_colors[system_size_index],
                                           linestyle="dashed", label=fr"$N = {length} \! \times \! {length}$")
            legend = axes[algorithm_index].legend(loc="lower left", fontsize=15.5, labelspacing=0)
            legend.get_frame().set_edgecolor("k")
            legend.get_frame().set_lw(2.5)

    observable_box_coords = [0.048, 0.55]
    axes[0].text(*observable_box_coords, r"$\Gamma = \epsilon^{-1}$ (electrolyte)", fontsize=14.5,
                 transform=axes[0].transAxes,
                 bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))
    axes[1].text(*observable_box_coords, r"$\Gamma = \Upsilon_\mathrm{HXY}$", fontsize=15.5,
                 transform=axes[1].transAxes,
                 bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))
    axes[2].text(*observable_box_coords, r"$\Gamma = \Upsilon_\mathrm{XY}$", fontsize=15.5,
                 transform=axes[2].transAxes,
                 bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))

    fig_label_x_coord = 1.46
    axes[0].text(fig_label_x_coord, 0.925, "(a)", fontsize=25.0)
    axes[1].text(fig_label_x_coord, 0.925, "(b)", fontsize=25.0)
    axes[2].text(fig_label_x_coord, 0.83, "(c)", fontsize=25.0)

    fig.savefig(f"{output_directory}/helicities_with_global_moves.pdf", bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


def make_topological_susceptibility_plots(algorithms, observables, model_temperatures, reduced_model_temperatures,
                                          linear_system_sizes, no_of_samples, no_of_runs, output_directory,
                                          sample_directories_by_algo, external_global_moves_string_all_moves,
                                          observable_plotting_markers, system_size_plotting_colors):
    fig, axes = make_three_empty_subfigures(r"$2 \pi^2 \beta J {\rm Var}[{\rm X}]$", 25.0)
    [axis.set_xlim(0.7, 1.6) for axis in axes]
    axes[1].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    axes[2].set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    for algorithm_index, algorithm in enumerate(algorithms):
        for system_size_index, length in enumerate(linear_system_sizes):
            no_of_sites_string = f"{length}x{length}_sites"
            """first compute mean acceptance rates across the runs"""
            get_means_and_errors("acceptance_rates", algorithm, model_temperatures[algorithm_index],
                                 no_of_sites_string, no_of_samples, external_global_moves_string_all_moves,
                                 output_directory, sample_directories_by_algo[algorithm_index][system_size_index],
                                 no_of_runs[algorithm_index])
            for observable_index, observable in enumerate(observables[algorithm_index]):
                if observable_index > 0:
                    means, errors = get_means_and_errors(observable, algorithm, model_temperatures[algorithm_index],
                                                         no_of_sites_string, no_of_samples,
                                                         external_global_moves_string_all_moves, output_directory,
                                                         sample_directories_by_algo[algorithm_index][system_size_index],
                                                         no_of_runs[algorithm_index])
                    if observable_index == 1:
                        axes[algorithm_index].errorbar(reduced_model_temperatures[algorithm_index], means, errors,
                                                       marker=observable_plotting_markers[observable_index - 1],
                                                       color=system_size_plotting_colors[system_size_index],
                                                       markersize=10, linestyle="dashed",
                                                       label=fr"$N = {length} \! \times \! {length}$")
                    elif observable_index == 2:
                        axes[algorithm_index].errorbar(reduced_model_temperatures[algorithm_index], means, errors,
                                                       marker=observable_plotting_markers[observable_index - 1],
                                                       color=system_size_plotting_colors[system_size_index],
                                                       markersize=10, linestyle="dashed")

        if algorithm_index == 0:
            legend = axes[algorithm_index].legend(loc="upper left", fontsize=15.0, labelspacing=0)
        else:
            legend = axes[algorithm_index].legend(handler_map={plt.Line2D: HandlerLine2D(update_func=remove_markers)},
                                                  loc="upper left", fontsize=15.0, labelspacing=0)
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_lw(2.5)

    axes[0].text(0.048, 0.44, r"$X = \mathbf{w}$ (electrolyte)", fontsize=14.0, transform=axes[0].transAxes,
                 bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))
    axes[1].text(0.048, 0.36, r"dots: $\,\, X = \mathbf{w}$ (HXY)" + "\n" +
                 r"stars: $X = \tilde{\mathbf{t}}^\mathrm{hot}$ (HXY)", fontsize=14.0, transform=axes[1].transAxes,
                 bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))
    axes[2].text(0.048, 0.36, r"dots: $\,\, X = \mathbf{w}$ (XY)" + "\n" +
                 r"stars: $X = \tilde{\mathbf{t}}^\mathrm{hot}$ (XY)", fontsize=14.0, transform=axes[2].transAxes,
                 bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))

    fig_label_x_coord = 1.49
    axes[0].text(fig_label_x_coord, 0.0, "(a)", fontsize=25.0)
    axes[1].text(fig_label_x_coord, 0.0, "(b)", fontsize=25.0)
    axes[2].text(fig_label_x_coord, 0.0, "(c)", fontsize=25.0)

    fig.savefig(f"{output_directory}/topological_susceptibilities_with_global_moves.pdf", bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


def make_topological_stability_plots(algorithms, observables, model_temperatures, reduced_model_temperatures,
                                     linear_system_sizes, no_of_samples, no_of_runs, output_directory,
                                     sample_directories_by_algo_local_moves, sample_directories_by_algo_all_moves,
                                     external_global_moves_string_local_moves, external_global_moves_string_all_moves,
                                     observable_plotting_markers, system_size_plotting_colors):
    fig, axes = make_three_empty_subfigures(
        r"$g_{\rm X}^{\rm local} \left( \tau = 10^6 \Delta t_\mathrm{Metrop}^\mathrm{model} \right)$",
        22.0, True)
    [axis.set_ylim(-0.49, 1.075) for axis in axes]
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

            """for our records, compute the means and errors of the inverse vacuum permittivity of each XY model"""
            if algorithm_index > 0:
                get_means_and_errors(
                    "inverse_vacuum_permittivity", algorithm, model_temperatures[algorithm_index],
                    no_of_sites_string, no_of_samples, external_global_moves_string_local_moves, output_directory,
                    sample_directories_by_algo_local_moves[algorithm_index][system_size_index],
                    no_of_runs[algorithm_index])
                get_means_and_errors(
                    "inverse_vacuum_permittivity", algorithm, model_temperatures[algorithm_index],
                    no_of_sites_string, no_of_samples, external_global_moves_string_all_moves, output_directory,
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
                if observable_index == 1:
                    axes[algorithm_index].plot(reduced_model_temperatures[algorithm_index], topological_stabilities,
                                               marker=observable_plotting_markers[observable_index - 1], markersize=10,
                                               color=system_size_plotting_colors[system_size_index], linestyle="dashed",
                                               label=fr"$N = {length} \! \times \! {length}$")
                elif observable_index == 2:
                    axes[algorithm_index].plot(reduced_model_temperatures[algorithm_index], topological_stabilities,
                                               marker=observable_plotting_markers[observable_index - 1], markersize=10,
                                               color=system_size_plotting_colors[system_size_index], linestyle="dashed")

        if algorithm_index == 0:
            legend = axes[algorithm_index].legend(loc="upper right", fontsize=15.0, labelspacing=0)
        else:
            legend = axes[algorithm_index].legend(handler_map={plt.Line2D: HandlerLine2D(update_func=remove_markers)},
                                                  loc="upper right", fontsize=15.0, labelspacing=0)
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_lw(2.5)

    text_box_coords = [0.04, 0.075]
    axes[0].text(*text_box_coords, r"$X = \mathbf{w}$ (electrolyte)", fontsize=14.0, transform=axes[0].transAxes,
                 bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))
    axes[1].text(*text_box_coords, r"dots: $\,\, X = \mathbf{w}$ (HXY)" + "\n" +
                 r"stars: $X = \tilde{\mathbf{t}}^\mathrm{hot}$ (HXY)", fontsize=14.0, transform=axes[1].transAxes,
                 bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))
    axes[2].text(*text_box_coords, r"dots: $\,\, X = \mathbf{w}$ (XY)" + "\n" +
                 r"stars: $X = \tilde{\mathbf{t}}^\mathrm{hot}$ (XY)", fontsize=14.0, transform=axes[2].transAxes,
                 bbox=dict(facecolor='none', edgecolor='black', linewidth=2.5, boxstyle='round, pad=0.5'))

    fig_labels = ["(a)", "(b)", "(c)"]
    fig_label_coords = [1.475, -0.4125]
    [axes[fig_index].text(*fig_label_coords, fig_labels[fig_index], fontsize=25.0) for fig_index in range(3)]

    fig.savefig(f"{output_directory}/topological_stabilities.pdf", bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


def make_three_empty_subfigures(y_axis_label, y_axis_label_font_size, remove_secondary_y_axis_values=False):
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))
    if remove_secondary_y_axis_values:
        axes[1].axes.yaxis.set_ticklabels([]), axes[2].axes.yaxis.set_ticklabels([])
        fig.tight_layout(w_pad=-0.5)
    else:
        fig.tight_layout(w_pad=1.4)
    axes[0].set_ylabel(y_axis_label, fontsize=y_axis_label_font_size, labelpad=-1.0)
    for axis_index, axis in enumerate(axes):
        axis.tick_params(which='both', direction='in', width=3)
        axis.tick_params(which='major', length=7, labelsize=25.0, pad=3)
        axis.tick_params(which='minor', length=4)
        [axis.spines[spine].set_linewidth(3.75) for spine in ["top", "bottom", "left", "right"]]
        if axis_index < 2:
            axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT}^{\rm elec} / \beta$", fontsize=25, labelpad=-0.5)
        else:
            axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta$", fontsize=25, labelpad=-0.5)
    return fig, axes


def remove_markers(handle, orig):
    handle.update_from(orig)
    handle.set_marker("")


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
