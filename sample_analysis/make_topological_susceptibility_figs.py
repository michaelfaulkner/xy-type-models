from make_summary_statistic_vs_temperature_plot import get_means_and_errors
import importlib
import matplotlib
import matplotlib.pyplot as plt
import os
import sample_getter
import sys

# import run script - have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def main():
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))
    fig.tight_layout(w_pad=2.5)
    for axis_index, axis in enumerate(axes):
        axis.tick_params(which='both', direction='in', width=3)
        axis.tick_params(which='major', length=7, labelsize=22.5, pad=3)
        axis.tick_params(which='minor', length=4)
        if axis_index < 2:
            axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT} / \beta$", fontsize=20, labelpad=-0.5)
        else:
            axis.set_xlabel(r"$\widetilde{\beta}_{\rm BKT}^{\rm XY} / \beta$", fontsize=20, labelpad=-0.5)
        '''axis.set_ylim(0.0, 1.0)
        axis.yaxis.set_major_locator(ticker.MultipleLocator(base=0.5))'''
        # axis.yaxis.set_major_formatter('{x:.1e}')
        if axis_index == 0:
            axis.set_ylabel(r"$\chi_{\rm w}$", fontsize=25, labelpad=-1.0)
        else:
            axis.set_ylabel(r"$\chi_{\rm x}$", fontsize=25, labelpad=-1.0)
            # axis.axes.yaxis.set_ticklabels([])
        [axis.spines[spine].set_linewidth(3.75) for spine in ["top", "bottom", "left", "right"]]

    algorithms = ["elementary-electrolyte", "hxy-uniform-noise-metropolis", "xy-uniform-noise-metropolis"]
    """Define observables[i] as the necessary observables of models[i]"""
    observables = [["topological_susceptibility"],
                   ["hxy_internal_twist_susceptibility", "potential_minimising_twist_susceptibility"],
                   ["xy_global_defect_susceptibility", "potential_minimising_twist_susceptibility"]]
    observable_plotting_markers = [".", "*"]
    get_sample_methods = [[getattr(sample_getter, "get_" + observable_string)
                           for observable_string in model_observables] for model_observables in observables]

    config_file_electrolyte = "config_files/top_susceptibilities/electrolyte_all.txt"
    config_file_hxy = "config_files/top_susceptibilities/hxy_all.txt"
    config_file_xy = "config_files/top_susceptibilities/xy_all.txt"
    (_, sample_directory_electrolyte, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps, no_of_observations,
     temperatures_electrolyte, use_external_global_moves, external_global_moves_string, no_of_runs, _, max_no_of_cpus
     ) = run_script.get_config_data(config_file_electrolyte)
    (_, sample_directory_hxy, _, _, _, _, temperatures_hxy, _, _, _, _, _) = run_script.get_config_data(config_file_hxy)
    (_, sample_directory_xy, _, _, _, _, temperatures_xy, _, _, _, _, _) = run_script.get_config_data(config_file_xy)
    output_directory = sample_directory_electrolyte.replace("/electrolyte_all", "")
    sample_directories = [sample_directory_electrolyte, sample_directory_hxy, sample_directory_xy]

    model_temperatures = [temperatures_electrolyte, temperatures_hxy, temperatures_xy]
    approx_transition_temperatures = [1.351, 1.351, 0.887]
    reduced_model_temperatures = [
        [temperature / approx_transition_temperatures[0] for temperature in temperatures_electrolyte],
        [temperature / approx_transition_temperatures[1] for temperature in temperatures_hxy],
        [temperature / approx_transition_temperatures[2] for temperature in temperatures_xy]]

    for algorithm_index, algorithm in enumerate(algorithms):
        for observable_index, observable in enumerate(observables[algorithm_index]):
            means, errors = get_means_and_errors(observable, algorithm, model_temperatures[algorithm_index],
                                                 no_of_sites, no_of_sites_string, no_of_equilibration_sweeps,
                                                 external_global_moves_string, sample_directories[algorithm_index],
                                                 no_of_runs, max_no_of_cpus)
            axes[algorithm_index].errorbar(reduced_model_temperatures[algorithm_index], means, errors,
                                           marker=observable_plotting_markers[observable_index], markersize=10,
                                           color="k")
        # plt.xlabel(r"reduced temperature, $\beta_{\rm c} / \beta$", fontsize=15, labelpad=10)
        # plt.ylabel(r"$\chi_{\rm{X}}$", fontsize=15, labelpad=10)

    fig.savefig(f"{output_directory}/top_susceptibilities_with_global_moves_{no_of_observations}_obs.pdf",
                bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: No positional argument allowed.")
    else:
        main()
