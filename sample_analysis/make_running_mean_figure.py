from sample_getter import get_physical_time_step
from setup_scripts import check_for_observable_error, get_sample_is_one_dimensional
import importlib
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sample_getter
import sys
import time

# import run script - have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def main(observable_string="rotated_magnetisation_phase"):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    config_file_2d = "config_files/running_mean_figure/2dxy.txt"
    config_file_3d = "config_files/running_mean_figure/3dxy.txt"
    (algorithm_name_2d, sample_directory_2d, no_of_sites_2d, no_of_sites_string_2d, no_of_equilibration_sweeps_2d, _,
     temperatures_2d, _, external_global_moves_string_2d, no_of_runs, _, _) = run_script.get_config_data(config_file_2d)
    (algorithm_name_3d, sample_directory_3d, no_of_sites_3d, no_of_sites_string_3d, no_of_equilibration_sweeps_3d, _,
     temperatures_3d, _, external_global_moves_string_3d, _, _, _) = run_script.get_config_data(config_file_3d)
    output_directory = sample_directory_2d.replace("/2dxy", "")

    check_for_observable_error(algorithm_name_2d, observable_string)
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
    linestyles = ["solid", "dotted", "dashed", "dashdot", (0, (1, 1)), (0, (5, 10)), (0, (5, 1)), (0, (3, 1, 1, 1))]
    sample_is_one_dimensional = get_sample_is_one_dimensional(observable_string)

    figure, axes = plt.subplots(2, 1, figsize=(7.5, 7.1))
    # axes[1].set_xlabel(r"$\tau$", fontsize=20, labelpad=0)
    axes[1].set_xlabel(r"$t_{\rm Metrop} / \Delta t_{\rm Metrop}$", fontsize=20, labelpad=0)
    axes[0].set_ylabel(r"$\bar{\phi}_m(\tau)$ / $s_{\phi_m}^2$", fontsize=20, labelpad=-10)
    axes[1].set_ylabel(r"$\bar{\phi}_m(\tau)$ / $s_{\phi_m}^2$", fontsize=20, labelpad=-10)
    [axis.tick_params(which='major', width=3, length=7, labelsize=18) for axis in axes]
    [axis.locator_params(axis='x', nbins=4) for axis in axes]
    [axis.locator_params(axis='y', nbins=8) for axis in axes]
    axes[0].xaxis.set_tick_params(labelbottom=False)
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes]
    colors = ["black", "red"]

    start_time = time.time()
    for temperature_index, temperature in enumerate(temperatures_2d):
        print(f"Temperature = {temperature:.4f}")
        try:
            physical_time_step_2d = np.loadtxt(
                f"{output_directory}/physical_time_step_{algorithm_name_2d.replace('-', '_')}_"
                f"{external_global_moves_string_2d}_{no_of_sites_string_2d}_temp_eq_{temperature:.4f}.csv", dtype=float,
                delimiter=",")
        except IOError:
            physical_time_step_2d = sum([
                get_physical_time_step(algorithm_name_2d, f"{sample_directory_2d}/run_{run_index}", temperature_index)
                for run_index in range(no_of_runs)]) / no_of_runs
            np.savetxt(f"{output_directory}/physical_time_step_{algorithm_name_2d.replace('-', '_')}_"
                       f"{external_global_moves_string_2d}_{no_of_sites_string_2d}_temp_eq_{temperature:.4f}.csv",
                       np.atleast_1d(physical_time_step_2d), delimiter=",")

        for run_index in range(no_of_runs):
            try:
                # first, attempt to open file containing the sample, running mean and running variance...
                sample_2d, running_mean_2d, running_variance_2d = np.load(
                    f"{output_directory}/{observable_string}_sample_and_running_mean_and_variance_"
                    f"{algorithm_name_2d.replace('-', '_')}_{external_global_moves_string_2d}_{no_of_sites_string_2d}_"
                    f"temp_eq_{temperature:.4f}_run_{run_index}.npy")
            except IOError:
                # ...then compute the sample, running mean and running variance if the file does not exist
                sample_2d = get_sample_method(f"{sample_directory_2d}/run_{run_index}", temperature, temperature_index,
                                              no_of_sites_2d, no_of_equilibration_sweeps_2d)
                if not sample_is_one_dimensional:
                    sample_2d = sample_2d.transpose()[0]
                running_mean_2d = np.array([np.mean(sample_2d[:index + 1]) for index in range(len(sample_2d))])
                running_variance_2d = np.array([np.var(sample_2d[:index + 1], ddof=1) for index in
                                                range(len(sample_2d))])
                np.save(f"{output_directory}/{observable_string}_sample_and_running_mean_and_variance_"
                        f"{algorithm_name_2d.replace('-', '_')}_{external_global_moves_string_2d}_"
                        f"{no_of_sites_string_2d}_temp_eq_{temperature:.4f}_run_{run_index}.npy",
                        np.array([sample_2d, running_mean_2d, running_variance_2d]))

            """Previous versions converted MC time to physical time, using the following commented-out code.  I then 
            realised that the only difference (in terms of symmetry breaking) between the 2D and 3D cases is the 
            absence of a singular limit (wrt N and then h).  But their equivalence (otherwise) means that we don't need 
            to compare the 2D and 3D sample means vs time, so I abandoned this figure and decided not to invest time in 
            finding the optimal temperatures to compares physical time across all fours phases that encompass the 2D 
            and 3D systems."""

            '''if temperature_index == 0 and run_index == 0:
                max_total_physical_time = physical_time_step_2d * (len(sample_2d) - 1)
            reduced_running_variance_2d = running_variance_2d[:int(max_total_physical_time / physical_time_step_2d) + 1]
            reduced_running_mean_2d = (running_mean_2d[:int(max_total_physical_time / physical_time_step_2d) + 1] /
                                       reduced_running_variance_2d[-1])
            if run_index == 0:
                axes[0].plot(np.arange(len(reduced_running_mean_2d)) * physical_time_step_2d, reduced_running_mean_2d,
                             color=colors[temperature_index], linewidth=2, linestyle=linestyles[run_index],
                             label=fr"$1 / (\beta J)$ = {temperature:.2f}")
            else:
                axes[0].plot(np.arange(len(reduced_running_mean_2d)) * physical_time_step_2d, reduced_running_mean_2d,
                             color=colors[temperature_index], linewidth=2, linestyle=linestyles[run_index])'''

            if run_index == 0:
                axes[0].plot(np.arange(len(running_mean_2d)), running_mean_2d / running_variance_2d[-1],
                             color=colors[temperature_index],  linewidth=2, linestyle=linestyles[run_index],
                             label=fr"$1 / (\beta J)$ = {temperature:.2f}")
            else:
                axes[0].plot(np.arange(len(running_mean_2d)), running_mean_2d / running_variance_2d[-1],
                             color=colors[temperature_index],  linewidth=2, linestyle=linestyles[run_index])

    for temperature_index, temperature in enumerate(temperatures_3d):
        observable_string = "rotated_magnetisation_phase"
        get_sample_method = getattr(sample_getter, "get_rotated_magnetisation_phase")
        print(f"Temperature = {temperature:.4f}")
        try:
            physical_time_step_3d = np.loadtxt(
                f"{output_directory}/physical_time_step_{algorithm_name_3d.replace('-', '_')}_"
                f"{external_global_moves_string_3d}_{no_of_sites_string_3d}_temp_eq_{temperature:.4f}.csv", dtype=float,
                delimiter=",")
        except IOError:
            physical_time_step_3d = sum([
                get_physical_time_step(algorithm_name_3d, f"{sample_directory_3d}/run_{run_index}", temperature_index)
                for run_index in range(no_of_runs)]) / no_of_runs
            np.savetxt(
                f"{output_directory}/physical_time_step_{algorithm_name_3d.replace('-', '_')}_"
                f"{external_global_moves_string_3d}_{no_of_sites_string_3d}_temp_eq_{temperature:.4f}.csv",
                np.atleast_1d(physical_time_step_3d), delimiter=",")

        for run_index in range(no_of_runs):
            try:
                # first, attempt to open file containing the sample, running mean and running variance...
                sample_3d, running_mean_3d, running_variance_3d = np.load(
                    f"{output_directory}/{observable_string}_sample_and_running_mean_and_variance_"
                    f"{algorithm_name_3d.replace('-', '_')}_{external_global_moves_string_3d}_"
                    f"{no_of_sites_string_3d}_temp_eq_{temperature:.4f}_run_{run_index}.npy")
            except IOError:
                # ...then compute the sample, running mean and running variance if the file does not exist
                sample_3d = get_sample_method(f"{sample_directory_3d}/run_{run_index}", temperature, temperature_index,
                                              no_of_sites_3d, no_of_equilibration_sweeps_3d)
                if not sample_is_one_dimensional:
                    sample_3d = sample_3d.transpose()[0]
                running_mean_3d = np.array([np.mean(sample_3d[:index + 1]) for index in range(len(sample_3d))])
                running_variance_3d = np.array([np.var(sample_3d[:index + 1], ddof=1) for index in
                                                range(len(sample_3d))])
                np.save(f"{output_directory}/{observable_string}_sample_and_running_mean_and_variance_"
                        f"{algorithm_name_3d.replace('-', '_')}_{external_global_moves_string_3d}_"
                        f"{no_of_sites_string_3d}_temp_eq_{temperature:.4f}_run_{run_index}.npy",
                        np.array([sample_3d, running_mean_3d, running_variance_3d]))

            if run_index == 0:
                axes[1].plot(np.arange(len(running_mean_3d)), running_mean_3d / running_variance_3d[-1],
                             color=colors[temperature_index], linewidth=2, linestyle=linestyles[run_index],
                             label=fr"$1 / (\beta J)$ = {temperature:.2f}")
            else:
                axes[1].plot(np.arange(len(running_mean_3d)), running_mean_3d / running_variance_3d[-1],
                             color=colors[temperature_index], linewidth=2, linestyle=linestyles[run_index])

    figure.tight_layout()
    legends = [axes[0].legend(loc="upper left", fontsize=12),
               axes[1].legend(loc="lower right", fontsize=12, title = "2D / 3D case in upper / lower panel")]
    [legend.get_frame().set_edgecolor("k") for legend in legends]
    [legend.get_frame().set_lw(3) for legend in legends]
    figure.savefig(f"{output_directory}/{observable_string}_running_means_vs_time_2D_vs_3D.pdf", bbox_inches="tight")
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("InterfaceError: No positional arguments required.  You may provide the string of the "
                        "observable whose running sample mean you wish to plot (default is "
                        "'rotated_magnetisation_phase').")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be the the string of the observable whose running sample "
              "mean you wish to plot (default is 'rotated_magnetisation_phase'.")
        main(sys.argv[1])
    else:
        main()
