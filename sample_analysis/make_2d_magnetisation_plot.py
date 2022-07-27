from markov_chain_diagnostics import get_cumulative_distribution
from sample_getter import get_cartesian_magnetisation, get_external_global_move, get_phase_in_polar_coordinates
from setup_scripts import check_initial_job_index, setup_pool
import importlib
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import sys
import time

# import run script - have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def main(config_file, no_of_histogram_bins=100):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps, no_of_observations,
     temperatures, use_external_global_moves, external_global_moves_string, no_of_jobs, initial_job_index,
     max_no_of_cpus) = run_script.get_config_data(config_file)
    if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
        print("ConfigurationError: The configuration file corresponds to a Maggs-electrolyte model but this script "
              "requires the XY of HXY model.")
        raise SystemExit
    check_initial_job_index(initial_job_index)

    start_time = time.time()
    if no_of_jobs == 1:
        make_plots(algorithm_name, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps,
                   temperatures, use_external_global_moves, external_global_moves_string, no_of_observations,
                   no_of_histogram_bins, job_index=None)
    else:
        pool = setup_pool(no_of_jobs, max_no_of_cpus)
        pool.starmap(make_plots, [
            (algorithm_name, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps,
             temperatures, use_external_global_moves, external_global_moves_string, no_of_observations,
             no_of_histogram_bins, job_index) for job_index in range(no_of_jobs)])
        pool.close()
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


def make_plots(algorithm_name, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps,
               temperatures, use_external_global_moves, external_global_moves_string, no_of_observations,
               no_of_histogram_bins, job_index):
    if job_index is None:
        job_index = 0
        sample_directory = output_directory
    else:
        sample_directory = f"{output_directory}/job_{job_index}"

    for temperature_index, temperature in enumerate(temperatures):
        print(f"Temperature = {temperature:.4f}")
        cartesian_magnetisation = get_cartesian_magnetisation(sample_directory, temperature, temperature_index,
                                                              no_of_sites)
        magnetisation_norm = np.linalg.norm(cartesian_magnetisation, axis=1)
        magnetisation_phase = np.array([get_phase_in_polar_coordinates(observation) for observation in
                                        cartesian_magnetisation])
        cartesian_magnetisation = cartesian_magnetisation.transpose()
        figure, axes = plt.subplots(1, 2, figsize=(20, 10))

        set_magnetisation_evolution_axes(axes)
        axes[0].plot(cartesian_magnetisation[0, :10000], cartesian_magnetisation[1, :10000], linestyle="solid",
                     color="black")
        axes[1].hist(magnetisation_phase[:10000], bins=no_of_histogram_bins, density=True, color="red", edgecolor="k")
        figure.savefig(f"{output_directory}/magnetisation_evolution_{algorithm_name.replace('-', '_')}_"
                       f"{external_global_moves_string}_{no_of_sites_string}_{no_of_observations}_obs_temp_eq_"
                       f"{temperature:.4f}_job_{job_index}_first_1e4_steps.pdf", bbox_inches="tight")
        [axis.cla() for axis in axes]

        cartesian_magnetisation = cartesian_magnetisation[:, no_of_equilibration_sweeps:]
        magnetisation_norm = magnetisation_norm[no_of_equilibration_sweeps:]
        magnetisation_phase = magnetisation_phase[no_of_equilibration_sweeps:]

        if use_external_global_moves and np.mean(magnetisation_norm) > (2.0 * no_of_sites) ** (-1.0 / 16.0):
            external_global_move = get_external_global_move(
                sample_directory, temperature, temperature_index, no_of_sites)[no_of_equilibration_sweeps:]
            if np.any(external_global_move != 0):
                global_move_event_times = np.where(np.linalg.norm(external_global_move, axis=1).astype(int) != 0)[0]
                time_window = 1000
                for index in global_move_event_times:
                    set_magnetisation_evolution_axes(axes)
                    """plot magnetisation at normalised time of global-move event +- time_window"""
                    axes[0].plot(
                        cartesian_magnetisation[0, max(0, index - time_window):min(len(magnetisation_norm) - 1,
                                                                                   index + time_window + 1)],
                        cartesian_magnetisation[1, max(0, index - time_window):min(len(magnetisation_norm) - 1,
                                                                                   index + time_window + 1)],
                        linestyle="solid", linewidth=1, color="black")
                    """plot magnetisation at time of global-move event with a red dot"""
                    axes[0].plot(cartesian_magnetisation[0, index:index + 1],
                                 cartesian_magnetisation[1, index:index + 1], "ro", markersize=2)
                    """estimate CDF of magnetisation phase at normalised time of global-move event +- time_window"""
                    axes[1].hist(magnetisation_phase[max(0, index - time_window):min(len(magnetisation_norm) - 1,
                                                                                     index + time_window + 1)],
                                 bins=no_of_histogram_bins, density=True, color="red", edgecolor="black")
                    figure.savefig(f"{output_directory}/magnetisation_evolution_{algorithm_name.replace('-', '_')}_"
                                   f"{external_global_moves_string}_{no_of_sites_string}_{no_of_observations}_obs_temp_"
                                   f"eq_{temperature:.4f}_job_{job_index}_around_global_twist_at_time_step_{index}.pdf",
                                   bbox_inches="tight")

                    [axis.cla() for axis in axes]

                    set_magnetisation_evolution_axes(axes)
                    """plot magnetisation up to time of global-move event"""
                    axes[0].plot(cartesian_magnetisation[0, 0:index + 1], cartesian_magnetisation[1, 0:index + 1],
                                 linestyle="solid", linewidth=1, color="black")
                    """plot magnetisation at time of global-move event with a red dot"""
                    axes[0].plot(cartesian_magnetisation[0, index:index + 1],
                                 cartesian_magnetisation[1, index:index + 1], "ro")
                    """estimate CDF of magnetisation phase up to time of global-move event"""
                    axes[1].hist(magnetisation_phase[0:index], bins=no_of_histogram_bins, density=True,
                                 color="red", edgecolor="black")
                    figure.savefig(f"{output_directory}/magnetisation_evolution_{algorithm_name.replace('-', '_')}_"
                                   f"{external_global_moves_string}_{no_of_sites_string}_{no_of_observations}_obs_temp_"
                                   f"eq_{temperature:.4f}_job_{job_index}_up_to_global_twist_at_time_step_{index}.pdf",
                                   bbox_inches="tight")
                    [axis.cla() for axis in axes]

        set_magnetisation_evolution_axes(axes)
        axes[1].tick_params(axis='y', colors='red')
        axes[1].set_ylabel(r"$\pi \left( \phi_m = x \right)$", fontsize=20, labelpad=4, color="red")
        axes[0].plot(cartesian_magnetisation[0, :10000], cartesian_magnetisation[1, :10000], linestyle="solid",
                     linewidth=1, color="black")
        axes[1].hist(magnetisation_phase, bins=no_of_histogram_bins, density=True, color="red", edgecolor="black")
        cdf_axis = axes[1].twinx()
        cdf_axis.set_ylim(0.0, 1.0)
        cdf_axis.tick_params(which='major', width=2, length=7, labelsize=18, pad=10)
        cdf_axis.set_ylabel(r"$F_{\phi_m,n}(x)$", fontsize=20, labelpad=-30)
        cdf_axis.plot(*get_cumulative_distribution(magnetisation_phase), color="black", linewidth=2, linestyle="-")
        cdf_axis.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
        cdf_axis.yaxis.set_major_formatter('{x:.1f}')
        figure.savefig(f"{output_directory}/magnetisation_evolution_{algorithm_name.replace('-', '_')}_"
                       f"{external_global_moves_string}_{no_of_sites_string}_{no_of_observations}_obs_temp_eq_"
                       f"{temperature:.4f}_job_{job_index}.pdf", bbox_inches="tight")
        figure.clear()

        plt.close()


def set_magnetisation_evolution_axes(axes):
    [axis.tick_params(which='both', width=2) for axis in axes]
    [axis.tick_params(which='major', length=7, labelsize=18, pad=10) for axis in axes]
    [axis.tick_params(which='minor', length=4) for axis in axes]

    axes[0].set_xlim([-1.0, 1.0])
    axes[0].set_ylim([-1.0, 1.0])

    # Move bottom x-axis and left y-axis to centre, passing through (0, 0)
    axes[0].spines['bottom'].set_position('center')
    axes[0].spines['left'].set_position('center')

    # Eliminate upper and right axes
    axes[0].spines['top'].set_color('none')
    axes[0].spines['right'].set_color('none')

    # Show ticks in the left and lower axes only
    axes[0].xaxis.set_ticks_position('bottom')
    axes[0].yaxis.set_ticks_position('left')

    axes[0].spines['left'].set_linewidth(2)
    axes[0].spines['bottom'].set_linewidth(2)

    minor_ticks = np.arange(-0.75, 1.0, 0.25)
    axes[0].set_xticks([-1.0, 1.0])
    axes[0].xaxis.set_major_formatter('{x:.1f}')
    axes[0].set_xticks(minor_ticks, minor=True)
    axes[0].set_yticks([-1.0, 1.0])
    axes[0].yaxis.set_major_formatter('{x:.1f}')
    axes[0].set_yticks(minor_ticks, minor=True)

    axes[0].set_xlabel(r"$m_x$", fontsize=20)
    axes[0].set_ylabel(r"$m_y$", fontsize=20, rotation="horizontal")
    axes[0].xaxis.set_label_coords(1.0, 0.55)
    axes[0].yaxis.set_label_coords(0.55, 0.9725)

    axes[1].set_xlim([-math.pi, math.pi])
    axes[1].set_xticks(np.arange(-math.pi, math.pi + 0.5 * math.pi / 2, step=(0.5 * math.pi)))
    axes[1].set_xticklabels([r"$-\pi$", r"$-\pi / 2$", r"$0$", r"$\pi / 2$", r"$\pi$"])
    axes[1].set_xlabel(r"$x$", fontsize=20, labelpad=8)
    axes[1].yaxis.set_major_formatter('{x:.2f}')
    axes[1].set_ylabel(r"$\pi \left( \phi_m = x \right)$", fontsize=20, labelpad=4)
    [axes[1].spines[spine].set_linewidth(2) for spine in ["top", "bottom", "left", "right"]]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise Exception("InterfaceError: One positional argument required - give the configuration-file location.")
    else:
        print("One positional argument provided.  It must be the location of the configuration file.")
        main(sys.argv[1])
