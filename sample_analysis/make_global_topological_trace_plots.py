from setup_scripts import check_initial_run_index, setup_pool
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
    (_, sample_directory_electrolyte, no_of_sites, no_of_sites_string,
     no_of_equilibration_sweeps, no_of_observations, temperatures_electrolyte, use_external_global_moves,
     external_global_moves_string, no_of_runs, initial_run_index, max_no_of_cpus
     ) = run_script.get_config_data(config_file_electrolyte)
    (_, sample_directory_hxy, _, _, _, _, temperatures_hxy, _, _, _, _, _) = run_script.get_config_data(config_file_hxy)
    (_, sample_directory_xy, _, _, _, _, temperatures_xy, _, _, _, _, _) = run_script.get_config_data(config_file_xy)
    output_directory = sample_directory_electrolyte.replace("/electrolyte", "")
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
        (sample_directories, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps,
         reduced_model_temperatures, no_of_observations, run_index) for run_index in range(no_of_runs)])
    pool.close()
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


def make_plots(sample_directories, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps,
               reduced_model_temperatures, no_of_observations, run_index):
    fig, axes = plt.subplots(2, 3, figsize=(30, 10))
    fig.tight_layout(w_pad=6.5)
    setup_figure_axes(axes)
    models = ["electrolyte", "hxy_model", "xy_model"]
    """Define observables[i] as the necessary observables of models[i]"""
    observables = [["electric_field_zero_mode", "topological_sector"],
                   ["macro_josephson_current", "hxy_topological_sector", "potential_minimising_twists"],
                   ["macro_josephson_current", "xy_twist_relaxation_field", "potential_minimising_twists",
                    "xy_topological_sector"]]
    observable_plotting_colours = ["black", "blue", "red", "yellow"]
    get_sample_methods = [[getattr(sample_getter, "get_" + observable_string)
                           for observable_string in model_observables] for model_observables in observables]
    run_indexed_sample_directories = [f"{sample_directory}/run_{run_index}" for sample_directory in sample_directories]

    for model_index, model in enumerate(models):
        for temperature_index, temperature in enumerate(reduced_model_temperatures[model_index]):
            samples = np.array([get_sample_method(run_indexed_sample_directories[model_index], temperature,
                                                  temperature_index, no_of_sites, no_of_equilibration_sweeps)
                                for get_sample_method in get_sample_methods[model_index]])
            if model_index < 2:
                samples[0] = (no_of_sites ** 0.5) * samples[0] / 2.0 / math.pi
            [axes[temperature_index, model_index].plot(samples[observable_index, :10000, 1], linestyle="solid",
                                                       color=observable_plotting_colours[observable_index])
             for observable_index in reversed(range(len(samples)))]
    fig.savefig(f"{output_directory}/global_topological_trace_plots_sans_global_moves_{no_of_sites_string}_"
                f"{no_of_observations}_obs_run_{run_index}.pdf", bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


def setup_figure_axes(axes):
    [axis.tick_params(which='both', direction='in', width=3) for axis in axes.flatten()]
    [axis.tick_params(which='major', length=7, labelsize=22.5, pad=3) for axis in axes.flatten()]
    [axis.tick_params(which='minor', length=4) for axis in axes.flatten()]
    [axis.spines[spine].set_linewidth(3.75) for spine in ["top", "bottom", "left", "right"] for axis in axes.flatten()]
    [axis.set_xlim([0.0, 5.0e3]) for axis in axes.flatten()]
    [axis.set_ylim([-1.5, 1.5]) for axis in axes.flatten()]
    [axes[1, model_index].set_xlabel(r"$t / \Delta t_{\mathrm{Metrop}}$", fontsize=25) for model_index in range(3)]
    for model_index in range(3):
        if model_index == 0:
            [axes[temp_index, model_index].set_ylabel(r"$\chi_{\rm w}(\widetilde{\beta}_{\rm BKT} / \beta = ?)$",
                                                      fontsize=25, labelpad=-1.0) for temp_index in range(2)]
        elif model_index == 1:
            [axes[temp_index, model_index].set_ylabel(r"$\chi_{\rm x}(\widetilde{\beta}_{\rm BKT} / \beta = ?)$",
                                                      fontsize=25, labelpad=-1.0) for temp_index in range(2)]
        else:
            [axes[temp_index, model_index].set_ylabel(
                r"$\chi_{\rm x}(\widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta = ?)$", fontsize=25, labelpad=-1.0)
                for temp_index in range(2)]


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: No positional argument allowed.")
    else:
        main()
