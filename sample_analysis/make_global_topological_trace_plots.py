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
     no_of_equilibration_sweeps, no_of_observations, temperatures, use_external_global_moves,
     external_global_moves_string, no_of_runs, initial_run_index, max_no_of_cpus
     ) = run_script.get_config_data(config_file_electrolyte)
    sample_directory_hxy = run_script.get_config_data(config_file_hxy)[1]
    sample_directory_xy = run_script.get_config_data(config_file_xy)[1]
    output_directory = sample_directory_electrolyte.replace("/electrolyte", "")
    sample_directories = [sample_directory_electrolyte, sample_directory_hxy, sample_directory_xy]

    check_initial_run_index(initial_run_index)
    start_time = time.time()
    pool = setup_pool(no_of_runs, max_no_of_cpus)
    pool.starmap(make_plots, [
        (sample_directories, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps,
         temperatures, no_of_observations, run_index) for run_index in range(no_of_runs)])
    pool.close()
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


def make_plots(sample_directories, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps,
               temperatures, no_of_observations, run_index):
    fig, axes = plt.subplots(2, 3, figsize=(20, 10))
    fig.tight_layout()
    setup_figure_axes(axes)
    models = ["electrolyte", "hxy_model", "xy_model"]
    """Define observables[i] as the necessary observables of models[i]"""
    observables = [["electric_field_zero_mode", "topological_sector"],
                   ["total_vortex_polarisation", "hxy_topological_sector", "potential_minimising_twists"],
                   ["total_vortex_polarisation", "xy_twist_relaxation_field", "potential_minimising_twists",
                    "xy_topological_sector"]]
    observable_plotting_colours = ["black", "blue", "red", "yellow"]
    get_sample_methods = [[getattr(sample_getter, "get_" + observable_string)
                           for observable_string in model_observables] for model_observables in observables]
    run_indexed_sample_directories = [f"{sample_directory}/run_{run_index}" for sample_directory in sample_directories]
    for temperature_index, temperature in enumerate(temperatures):
        samples = [[get_sample_method(run_indexed_sample_directories[model_index], temperature, temperature_index,
                                      no_of_sites, no_of_equilibration_sweeps)
                    for get_sample_method in get_sample_methods[model_index]] for model_index in range(3)]
        [[axes[temperature_index, model_index].plot(
            list(zip(*samples[model_index][observable_index][:10000]))[1], linestyle="solid",
            color=observable_plotting_colours[observable_index]) for observable_index in range(len(model_observables))]
            for model_index, model_observables in enumerate(observables)]
    fig.savefig(f"{output_directory}/global_topological_trace_plots_sans_global_moves_{no_of_sites_string}_"
                f"{no_of_observations}_obs_run_{run_index}.pdf", bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


def setup_figure_axes(axes):
    [axis.tick_params(which='both', width=2) for axis in axes.flatten()]
    [axis.tick_params(which='major', length=7, labelsize=18, pad=10) for axis in axes.flatten()]
    [axis.tick_params(which='minor', length=4) for axis in axes.flatten()]
    # [axis.set_linewidth(2) for axis in axes.flatten()]

    [axes[1, model_index].set_xlim([0.0, 5.0e3]) for model_index in range(3)]
    [axis.set_ylim([-1.5, 1.5]) for axis in axes.flatten()]
    [axes[1, model_index].set_xlabel(r"$t$", fontsize=20) for model_index in range(3)]
    [axis.set_ylabel(r"???", fontsize=20, rotation="horizontal") for axis in axes.flatten()]

    # todo change the following ticks parameters
    """[axes[1, model_index].set_xticks(np.arange(-math.pi, math.pi + 0.5 * math.pi / 2, step=(0.5 * math.pi)))
     for model_index in range(3)]
    [axes[1, model_index].set_xticklabels([r"set", r"these"]) for model_index in range(3)]"""


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: No positional argument allowed.")
    else:
        main()
