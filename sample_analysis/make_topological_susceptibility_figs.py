from make_summary_statistic_vs_temperature_plot import get_means_and_errors
from setup_scripts import setup_pool
import importlib
import matplotlib
import matplotlib.pyplot as plt
import os
import sys

# import run script - have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def main(no_of_system_sizes=5):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"

    linear_system_sizes = [2 ** (index + 2) for index in range(no_of_system_sizes)]
    base_config_file_electrolyte_all_moves = f"config_files/top_susc_figs/4x4_elec_all.txt"
    base_config_file_hxy_all_moves = f"config_files/top_susc_figs/4x4_hxy_all.txt"
    base_config_file_xy_all_moves = f"config_files/top_susc_figs/4x4_xy_all.txt"
    (algorithm_name_electrolyte, sample_directory_4x4_electrolyte_all_moves, _, _, no_of_equilibration_sweeps,
     no_of_observations, temperatures_electrolyte, _, external_global_moves_string_all_moves,
     no_of_runs, _, max_no_of_cpus) = run_script.get_config_data(base_config_file_electrolyte_all_moves)
    (algorithm_name_hxy, _, _, _, _, _, temperatures_hxy, _, _, _, _, _
     ) = run_script.get_config_data(base_config_file_hxy_all_moves)
    (algorithm_name_xy, _, _, _, _, _, temperatures_xy, _, _, _, _, _
     ) = run_script.get_config_data(base_config_file_xy_all_moves)

    algorithms = [algorithm_name_electrolyte, algorithm_name_hxy, algorithm_name_xy]
    """Define observables[i] as the necessary observables of algorithms[i]"""
    observables = [["topological_susceptibility"],
                   ["hxy_internal_twist_susceptibility", "potential_minimising_twist_susceptibility"],
                   ["xy_global_defect_susceptibility", "potential_minimising_twist_susceptibility"]]
    observable_plotting_markers = [".", "*"]
    system_size_plotting_colors = ["black", "red", "blue", "green", "magenta", "indigo"][:no_of_system_sizes]
    system_size_plotting_colors.reverse()

    output_directory = sample_directory_4x4_electrolyte_all_moves.replace("/4x4_elec_all", "")
    sample_directories_electrolyte_all_moves = [f"{output_directory}/{length}x{length}_elec_all" for length in
                                                linear_system_sizes]
    sample_directories_hxy_all_moves = [f"{output_directory}/{length}x{length}_hxy_all" for length in
                                        linear_system_sizes]
    sample_directories_xy_all_moves = [f"{output_directory}/{length}x{length}_xy_all" for length in linear_system_sizes]
    sample_directories_by_algo = [sample_directories_electrolyte_all_moves, sample_directories_hxy_all_moves,
                                  sample_directories_xy_all_moves]
    model_temperatures = [temperatures_electrolyte, temperatures_hxy, temperatures_xy]
    approx_transition_temperatures = [1.351, 1.351, 0.887]
    reduced_model_temperatures = [
        [temperature / approx_transition_temperatures[0] for temperature in temperatures_electrolyte],
        [temperature / approx_transition_temperatures[1] for temperature in temperatures_hxy],
        [temperature / approx_transition_temperatures[2] for temperature in temperatures_xy]]

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

    for algorithm_index, algorithm in enumerate(algorithms):
        for system_size_index, length in enumerate(linear_system_sizes):
            no_of_sites_string = f"{length}x{length}_sites"
            """first compute mean acceptance rates across the runs"""
            get_means_and_errors("acceptance_rates", algorithm, model_temperatures[algorithm_index], length ** 2,
                                 no_of_sites_string, no_of_equilibration_sweeps,
                                 external_global_moves_string_all_moves,
                                 sample_directories_by_algo[algorithm_index][system_size_index], no_of_runs,
                                 max_no_of_cpus)
            for observable_index, observable in enumerate(observables[algorithm_index]):
                means, errors = get_means_and_errors(observable, algorithm, model_temperatures[algorithm_index],
                                                     length ** 2, no_of_sites_string, no_of_equilibration_sweeps,
                                                     external_global_moves_string_all_moves,
                                                     sample_directories_by_algo[algorithm_index][system_size_index],
                                                     no_of_runs, max_no_of_cpus)
                axes[algorithm_index].errorbar(reduced_model_temperatures[algorithm_index], means, errors,
                                               marker=observable_plotting_markers[observable_index], markersize=10,
                                               color=system_size_plotting_colors[system_size_index])
    fig.savefig(f"{output_directory}/top_susceptibilities_with_global_moves_{no_of_observations}_obs.pdf",
                bbox_inches="tight")
    [axis.cla() for axis in axes.flatten()]
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("InterfaceError: At most one positional arguments permitted.  None are required but you may "
                        "provide no_of_system_sizes, which must be an integer greater than 0 and less than 6 (default "
                        "value is 5).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be no_of_system_sizes - which must be an integer greater "
              "than 0 and less than 6 (default value is 5).")
        chosen_no_of_system_sizes = int(sys.argv[1])
        if chosen_no_of_system_sizes < 1 or chosen_no_of_system_sizes > 5:
            raise Exception("InterfaceError: no_of_system_sizes must be an integer greater than 0 and less than 6 "
                            "(default value is 5).")
        main(chosen_no_of_system_sizes)
    else:
        print("No positional arguments provided.  None are required but you may provide no_of_system_sizes, which must "
              "be an integer greater than 0 and less than 6 (default value is 5).")
        main()
