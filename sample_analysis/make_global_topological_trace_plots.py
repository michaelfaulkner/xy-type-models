from setup_scripts import check_initial_run_index, reverse_enumerate, setup_pool
from sample_getter import get_macro_josephson_current
import sample_getter
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


def main():
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    config_file_electrolyte = "config_files/global_top_trace_plots/electrolyte.txt"
    config_file_hxy = "config_files/global_top_trace_plots/hxy.txt"
    config_file_xy = "config_files/global_top_trace_plots/xy.txt"
    (algorithm_name_electrolyte, sample_directory_electrolyte, no_of_sites, no_of_sites_string,
     no_of_equilibration_sweeps, no_of_observations, temperatures_electrolyte, use_external_global_moves,
     external_global_moves_string, no_of_runs, initial_run_index, max_no_of_cpus
     ) = run_script.get_config_data(config_file_electrolyte)
    (algorithm_name_hxy, sample_directory_hxy, _, _, _, _, temperatures_hxy, _, _, _, _, _
     ) = run_script.get_config_data(config_file_hxy)
    (algorithm_name_xy, sample_directory_xy, _, _, _, _, temperatures_xy, _, _, _, _, _
     ) = run_script.get_config_data(config_file_xy)
    output_directory = sample_directory_electrolyte.replace("/electrolyte", "")
    algorithm_names = [algorithm_name_electrolyte, algorithm_name_hxy, algorithm_name_xy]
    sample_directories = [sample_directory_electrolyte, sample_directory_hxy, sample_directory_xy]
    approx_transition_temperatures = [1.351, 1.351, 0.887]
    reduced_model_temperatures = [
        [temperature / approx_transition_temperatures[0] for temperature in temperatures_electrolyte],
        [temperature / approx_transition_temperatures[1] for temperature in temperatures_hxy],
        [temperature / approx_transition_temperatures[2] for temperature in temperatures_xy]]

    check_initial_run_index(initial_run_index)
    start_time = time.time()
    pool = setup_pool(no_of_runs, max_no_of_cpus)
    pool.starmap(make_plots, [
        (algorithm_names, sample_directories, output_directory, no_of_sites, no_of_sites_string,
         reduced_model_temperatures, no_of_equilibration_sweeps, no_of_observations, run_index)
        for run_index in range(no_of_runs)])
    pool.close()
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


def make_plots(algorithm_names, sample_directories, output_directory, no_of_sites, no_of_sites_string,
               reduced_model_temperatures, no_of_equilibration_sweeps, no_of_observations, run_index):
    fig, axes = plt.subplots(2, 3, figsize=(30, 10))
    fig.tight_layout(w_pad=4.25)
    left_letter_width = 0.299
    middle_letter_width = 0.635
    right_letter_width = 0.9725
    letter_heights = 0.94
    fig.text(left_letter_width, letter_heights, "(a)", fontsize=27)
    fig.text(middle_letter_width, letter_heights, "(b)", fontsize=27)
    fig.text(right_letter_width, letter_heights, "(c)", fontsize=27)
    letter_heights = 0.445
    fig.text(left_letter_width, letter_heights, "(d)", fontsize=27)
    fig.text(middle_letter_width, letter_heights, "(e)", fontsize=27)
    fig.text(right_letter_width, letter_heights, "(f)", fontsize=27)
    setup_figure_axes(axes)
    models = ["electrolyte", "hxy_model", "xy_model"]
    """Define observables[i] as the observables of interest of models[i]..."""
    observables = [["electric_field_zero_mode", "topological_sector"],
                   ["macro_josephson_current", "hxy_topological_sector", "potential_minimising_twists"],
                   ["xy_emergent_field_zero_mode", "xy_topological_sector", "potential_minimising_twists",
                    "xy_twist_relaxation_field"]]
    """...then define the Cartesian coordinate of each observable (x = 0 / y = 1 for electro/spin-rep quantities)..."""
    observable_cartesian_coords = [[0, 0], [1, 0, 1], [0, 0, 1, 1]]
    """...then define the corresponding plot labels"""
    plot_labels = [[r"$\bar{E}_{\rm x} / (2\pi J / L)$", r"$w_{\rm x}$"],
                   [r"$\bar{E}_{\rm x} / (2\pi J / L)$", r"$w_{\rm x}$",
                    r"$\tilde{t}_{\rm y}^{\rm non-anneal}$"],
                   [r"$\bar{E}_{\rm x} / (2\pi J / L)$", r"$w_{\rm x}$", r"$\tilde{t}_{\rm y}^{\rm non-anneal}$",
                    r"$\tilde{t}_{\rm y}$"]]
    observable_plotting_colours = ["black", "blue", "red", "orange"]
    observable_marker_styles = [None, ".", "*", "D"]
    get_sample_methods = [[getattr(sample_getter, "get_" + observable_string)
                           for observable_string in model_observables] for model_observables in observables]
    run_indexed_sample_directories = [f"{sample_directory}/run_{run_index}" for sample_directory in sample_directories]

    for model_index, model in enumerate(models):
        for temp_index, temperature in enumerate(reduced_model_temperatures[model_index]):
            compressed_sample_filenames = [
                (f"{output_directory}/{observable_string}_{algorithm_names[model_index].replace('-', '_')}_"
                 f"sans_global_moves_{no_of_sites_string}_{no_of_observations}_obs_temp_eq_{temperature:.4f}_"
                 f"run_{run_index}.npy") for observable_string in observables[model_index]]
            try:
                samples = np.array([np.load(compressed_sample_filenames[observable_index])
                                    for observable_index in range(len(observables[model_index]))])
            except IOError:
                samples = np.array([get_sample_method(run_indexed_sample_directories[model_index], temperature,
                                                      temp_index, no_of_sites, no_of_equilibration_sweeps)
                                    for get_sample_method in get_sample_methods[model_index]])
                [np.save(compressed_sample_filenames[observable_index], samples[observable_index])
                 for observable_index in range(len(observables[model_index]))]
            samples[0] = (no_of_sites ** 0.5) * samples[0] / 2.0 / math.pi
            if temp_index == 0:
                for observable_index in reversed(range(len(samples))):
                    if observable_index == 0:
                        axes[temp_index, model_index].plot(samples[observable_index, :no_of_observations, 
                                                           observable_cartesian_coords[model_index][observable_index]],
                                                           linestyle="-", linewidth=2.5, 
                                                           color=observable_plotting_colours[observable_index],
                                                           label=f"X = {plot_labels[model_index][observable_index]}")
                    else:
                        axes[temp_index, model_index].plot(samples[observable_index, :no_of_observations,
                                                           observable_cartesian_coords[model_index][observable_index]],
                                                           marker=observable_marker_styles[observable_index],
                                                           markersize=12, linestyle="-",  
                                                           color=observable_plotting_colours[observable_index],
                                                           label=f"X = {plot_labels[model_index][observable_index]}")
            else:
                for observable_index in reversed(range(len(samples))):
                    if observable_index == 0:
                        axes[temp_index, model_index].plot(samples[observable_index, :no_of_observations,
                                                           observable_cartesian_coords[model_index][observable_index]],
                                                           linestyle="-", linewidth=2.5,
                                                           color=observable_plotting_colours[observable_index])
                    else:
                        axes[temp_index, model_index].plot(samples[observable_index, :no_of_observations,
                                                           observable_cartesian_coords[model_index][observable_index]],
                                                           marker=observable_marker_styles[observable_index],
                                                           markersize=12, linestyle="-",
                                                           color=observable_plotting_colours[observable_index])
    for model_index in range(len(models)):
        handles, labels = axes[0, model_index].get_legend_handles_labels()
        if model_index == 0:
            legend = axes[0, model_index].legend(reversed(handles), reversed(labels), loc="upper left", fontsize=20)
        else:
            legend = axes[0, model_index].legend(reversed(handles), reversed(labels), loc="upper left", fontsize=20,
                                                 ncol=2)
        legend.get_frame().set_edgecolor("k")
        legend.get_frame().set_lw(3)

    inset_axis = plt.axes((0.8375, 0.545, 0.15, 0.15))
    inset_axis.tick_params(which='both', direction='in', length=5, width=3, labelsize=12, pad=2)
    inset_axis.set_xlim([0.0, 2.499e3])
    inset_axis.xaxis.set_label_position("top")
    inset_axis.xaxis.tick_top()
    inset_axis.set_xlabel(r"$t / \Delta t_{\mathrm{Metrop}}$", fontsize=14, labelpad=3)
    inset_axis.set_ylim(-0.6, 0.6)
    inset_axis.set_ylabel(r"$\sum_{\langle \mathbf{r}, \mathbf{r}' \rangle_y} "
                          r"\sin (\varphi_\mathbf{r} - \varphi_{\mathbf{r}'}) / N$", fontsize=12, labelpad=-8)
    [inset_axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]
    for temp_index, temperature in reverse_enumerate(reduced_model_temperatures[2]):
        xy_macro_josephson_current_sample = get_macro_josephson_current(
            run_indexed_sample_directories[2], temperature, temp_index, no_of_sites, no_of_equilibration_sweeps)
        if temp_index == 0:
            inset_axis.plot(xy_macro_josephson_current_sample[:no_of_observations, 1], linestyle="solid",
                            color="black", label=r"$\widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta = 0.95$")
        else:
            inset_axis.plot(xy_macro_josephson_current_sample[:no_of_observations, 1], linestyle="solid",
                            color="red", label=r"$\widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta = 1.5$")
    handles, labels = inset_axis.get_legend_handles_labels()
    legend = inset_axis.legend(reversed(handles), reversed(labels), loc="lower left", fontsize=8, ncol=2)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(3)

    fig.savefig(f"{output_directory}/global_topological_trace_plots_sans_global_moves_{no_of_sites_string}_"
                f"{no_of_observations}_obs_run_{run_index}.pdf", bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


def setup_figure_axes(axes):
    [axis.tick_params(which='both', direction='in', width=3) for axis in axes.flatten()]
    [axis.tick_params(which='major', length=7, labelsize=22.5, pad=3) for axis in axes.flatten()]
    [axis.tick_params(which='minor', length=4) for axis in axes.flatten()]
    [axis.spines[spine].set_linewidth(3.75) for spine in ["top", "bottom", "left", "right"] for axis in axes.flatten()]
    [axis.set_xlim([0.0, 2.5e3]) for axis in axes.flatten()]
    [axis.set_ylim([-2.1, 2.1]) for axis in axes.flatten()]
    [axes[1, model_index].set_xlabel(r"$t / \Delta t_{\mathrm{Metrop}}$", fontsize=25) for model_index in range(3)]
    for model_index in range(3):
        if model_index == 0:
            for temp_index in range(2):
                if temp_index == 0:
                    axes[temp_index, model_index].set_ylabel(r"$X(t; \widetilde{\beta}_{\rm BKT} / \beta = 0.95)$",
                                                             fontsize=25, labelpad=-1.0)
                else:
                    axes[temp_index, model_index].set_ylabel(r"$X(t; \widetilde{\beta}_{\rm BKT} / \beta = 1.5)$",
                                                             fontsize=25, labelpad=-1.0)
        elif model_index == 1:
            for temp_index in range(2):
                if temp_index == 0:
                    axes[temp_index, model_index].set_ylabel(r"$X(t; \widetilde{\beta}_{\rm BKT} / \beta = 0.95)$",
                                                             fontsize=25, labelpad=-1.0)
                else:
                    axes[temp_index, model_index].set_ylabel(r"$X(t; \widetilde{\beta}_{\rm BKT} / \beta = 1.5)$",
                                                             fontsize=25, labelpad=-1.0)
        else:
            for temp_index in range(2):
                if temp_index == 0:
                    axes[temp_index, model_index].set_ylabel(
                        r"$X(t; \widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta = 0.95)$", fontsize=25, labelpad=-1.0)
                else:
                    axes[temp_index, model_index].set_ylabel(
                        r"$X(t; \widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta = 1.5)$", fontsize=25, labelpad=-1.0)


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: No positional argument allowed.")
    else:
        main()
