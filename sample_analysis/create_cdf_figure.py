import importlib
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import sys
import time

# import additional modules; have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")
setup_scripts = importlib.import_module("setup_scripts")
sample_getter = importlib.import_module("sample_getter")
run_script = importlib.import_module("run")


def main():
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    matplotlib.rcParams["axes.unicode_minus"] = False
    config_file_64x64_metrop = "config_files/cdf_figure/64x64_metropolis.txt"
    config_file_64x64_ecmc = "config_files/cdf_figure/64x64_ecmc.txt"
    config_file_256x256_metrop = "config_files/cdf_figure/256x256_metropolis.txt"
    config_file_256x256_ecmc = "config_files/cdf_figure/256x256_ecmc.txt"
    (algorithm_name_metrop, sample_directory_64x64_metrop, no_of_sites_64x64, no_of_sites_string_64x64,
     no_of_equilibration_sweeps_metrop, no_of_observations_metrop, temperatures, use_external_global_moves,
     external_global_moves_string, no_of_jobs_metrop, _, max_no_of_cpus) = run_script.get_config_data(
        config_file_64x64_metrop)
    (algorithm_name_ecmc, sample_directory_64x64_ecmc, _, _, no_of_equilibration_sweeps_ecmc, no_of_observations_ecmc,
     _, _, _, no_of_jobs_ecmc, _, _) = run_script.get_config_data(config_file_64x64_ecmc)
    (_, sample_directory_256x256_metrop, no_of_sites_256x256, no_of_sites_string_256x256, _, _, _, _, _, _, _, _
     ) = run_script.get_config_data(config_file_256x256_metrop)
    sample_directory_256x256_ecmc = run_script.get_config_data(config_file_256x256_ecmc)[1]
    output_directory = sample_directory_64x64_metrop.replace("/64x64_metropolis", "")

    figure, axes = plt.subplots(1, 3, figsize=(15, 4.25))
    alphabetic_labels = ["(a)", "(b)", "(c)"]
    [axis.text(-0.25, -0.25, f"{alphabetic_labels[axis_index]}", fontsize=18) for axis_index, axis in enumerate(axes)]
    # plt.subplots_adjust(wspace=0, hspace=0)
    figure.tight_layout(w_pad=-6.0)
    for axis_index, axis in enumerate(axes):
        axis.tick_params(which='both', direction='in', width=3)
        axis.tick_params(which='major', length=7, labelsize=18, pad=5)
        axis.tick_params(which='minor', length=4)
        axis.set_xlim([-math.pi, math.pi])
        axis.set_xticks(np.arange(-math.pi, math.pi + 0.5 * math.pi / 2, step=(0.5 * math.pi)))
        axis.set_xticklabels(["", r"$-\pi / 2$", r"$0$", r"$\pi / 2$", ""])
        axis.set_xlabel(r"$x$", fontsize=20, labelpad=-1)
        axis.set_ylim(0.0, 1.0)
        axis.yaxis.set_major_locator(ticker.MultipleLocator(base=0.5))
        if axis_index == 0:
            # axis.tick_params(which='major', width=3, length=7, labelsize=18, pad=10)
            axis.set_ylabel(r"$F_{\phi_m, n}(x)$", fontsize=20, labelpad=5)
            axis.yaxis.set_major_formatter('{x:.1f}')
        else:
            axis.axes.yaxis.set_ticklabels([])
        [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]

    inset_axis = plt.axes([0.0745, 0.71, 0.09, 0.22])
    inset_axis.tick_params(which='both', direction='in', length=5, width=3, labelsize=11)
    inset_axis.set_xlim([-math.pi, math.pi])
    inset_axis.set_xticks(np.arange(-math.pi, math.pi + 0.5 * math.pi / 2, step=(0.5 * math.pi)))
    inset_axis.set_xticklabels(["", r"$-\pi / 2$", r"$0$", r"$\pi / 2$", ""])
    inset_axis.set_xlabel(r"$x$", fontsize=12, labelpad=-1)
    inset_axis.set_ylim(0.0, 1.0)
    inset_axis.yaxis.set_label_position("right")
    inset_axis.yaxis.tick_right()
    inset_axis.set_ylabel(r"$F_{\phi_m, n}(x)$", fontsize=12, labelpad=-10)
    inset_axis.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    inset_axis.yaxis.set_major_formatter('{x:.1f}')
    inset_axis.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.5))
    [inset_axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]

    colors = ["red", "black"]
    linestyles = ["solid", "dotted", "dashed", "dashdot", (0, (1, 1)), (0, (5, 10)), (0, (5, 1)), (0, (3, 1, 1, 1))]

    start_time = time.time()
    make_non_schematic_subplot(algorithm_name_ecmc, algorithm_name_metrop, external_global_moves_string,
                               no_of_sites_64x64, no_of_sites_string_64x64, temperatures,
                               no_of_equilibration_sweeps_ecmc, no_of_equilibration_sweeps_metrop,
                               no_of_observations_ecmc, no_of_observations_metrop, no_of_jobs_ecmc, no_of_jobs_metrop,
                               output_directory, sample_directory_64x64_ecmc, sample_directory_64x64_metrop, axes,
                               inset_axis, linestyles, colors, 0)
    make_non_schematic_subplot(algorithm_name_ecmc, algorithm_name_metrop, external_global_moves_string,
                               no_of_sites_256x256, no_of_sites_string_256x256, temperatures,
                               no_of_equilibration_sweeps_ecmc, no_of_equilibration_sweeps_metrop,
                               no_of_observations_ecmc, no_of_observations_metrop, no_of_jobs_ecmc, no_of_jobs_metrop,
                               output_directory, sample_directory_256x256_ecmc, sample_directory_256x256_metrop, axes,
                               inset_axis, linestyles, colors, 1)
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    axes[2].plot([-math.pi, math.pi], [0.0, 1.0], linestyle="-", color="black", linewidth=2, label="symmetric phase")
    try:
        heaviside_centres = np.load(f"{output_directory}/heaviside_centres.npy")
    except IOError:
        heaviside_centres = [np.random.uniform(-math.pi, math.pi) for _ in range(8)]
        np.save(f"{output_directory}/heaviside_centres.npy", np.array(heaviside_centres))
    [axes[2].plot([heaviside_centre, heaviside_centre], [0.0, 1.0], color="red",
                  linestyle=linestyles[realisation_index], linewidth=2, label="broken-symmetry phase")
     for realisation_index, heaviside_centre in enumerate(heaviside_centres)]

    handles, labels = axes[0].get_legend_handles_labels()
    legends = [axes[0].legend(reversed(handles), reversed(labels), loc="lower right", fontsize=14),
               axes[1].legend(reversed(handles), reversed(labels), loc="upper left", fontsize=11)]
    [legend.get_frame().set_edgecolor("k") for legend in legends]
    [legend.get_frame().set_lw(3) for legend in legends]
    figure.savefig(f"{output_directory}/magnetisation_phase_cdfs.pdf", bbox_inches="tight")


def make_non_schematic_subplot(algorithm_name_ecmc, algorithm_name_metrop, external_global_moves_string, no_of_sites,
                               no_of_sites_string, temperatures, no_of_equilibration_sweeps_ecmc,
                               no_of_equilibration_sweeps_metrop, no_of_observations_ecmc, no_of_observations_metrop,
                               no_of_jobs_ecmc, no_of_jobs_metrop, output_directory, sample_directory_ecmc,
                               sample_directory_metrop, axes, inset_axis, linestyles, colors, subplot_index):
    try:
        for temperature_index, temperature in setup_scripts.reverse_enumerate(temperatures):
            for job_index in range(no_of_jobs_metrop):
                output_file_metrop = (
                    f"{output_directory}/magnetisation_phase_cdf_{algorithm_name_metrop.replace('-', '_')}_"
                    f"{external_global_moves_string}_{no_of_sites_string}_{no_of_observations_metrop}_obs_"
                    f"temp_eq_{temperature:.4f}_job_{job_index}.npy")
                if job_index == 0:
                    axes[subplot_index].plot(*np.load(output_file_metrop), color=colors[temperature_index],
                                             linestyle=linestyles[job_index], linewidth=2,
                                             label=fr"$1 / (\beta J)$ = {temperature:.2f}")
                    """use alternative CDF calculation in 
                        markov_chain_diagnostics.get_cumulative_distribution() to use fill_between() in the 
                        following line"""
                    """if temperature_index != 0:
                        axis.fill_between(cdf[0], cdf[1], np.arange(0.0, 1.0, 1.0 / float(len(cdf[0]))), 
                                          color='green')"""
                else:
                    axes[subplot_index].plot(*np.load(output_file_metrop), color=colors[temperature_index],
                                             linestyle=linestyles[job_index], linewidth=2)
                if job_index < no_of_jobs_ecmc:
                    output_file_ecmc = (
                        f"{output_directory}/magnetisation_phase_cdf_{algorithm_name_ecmc.replace('-', '_')}_"
                        f"{external_global_moves_string}_{no_of_sites_string}_{no_of_observations_ecmc}_obs_"
                        f"temp_eq_{temperature:.4f}_job_{job_index}.npy")
                    if subplot_index == 1:
                        inset_axis.plot(*np.load(output_file_ecmc), color=colors[temperature_index],
                                        linestyle=linestyles[job_index], linewidth=2)
    except IOError:
        for temperature_index, temperature in setup_scripts.reverse_enumerate(temperatures):
            print(f"Temperature = {temperature:.4f}")
            cdfs_of_magnetisation_phase_metrop = [get_cdf_of_magnetisation_phase(
                f"{sample_directory_metrop}/job_{job_index}", temperature, temperature_index, no_of_sites,
                no_of_equilibration_sweeps_metrop) for job_index in range(no_of_jobs_metrop)]
            cdfs_of_magnetisation_phase_ecmc = [get_cdf_of_magnetisation_phase(
                f"{sample_directory_ecmc}/job_{job_index}", temperature, temperature_index, no_of_sites,
                no_of_equilibration_sweeps_ecmc) for job_index in range(no_of_jobs_ecmc)]
            for job_index in range(no_of_jobs_metrop):
                if job_index == 0:
                    axes[subplot_index].plot(*cdfs_of_magnetisation_phase_metrop[job_index],
                                             color=colors[temperature_index], linestyle=linestyles[job_index],
                                             linewidth=2, label=fr"$1 / (\beta J)$ = {temperature:.2f}")
                else:
                    axes[subplot_index].plot(*cdfs_of_magnetisation_phase_metrop[job_index],
                                             color=colors[temperature_index], linestyle=linestyles[job_index],
                                             linewidth=2)
                np.save(f"{output_directory}/magnetisation_phase_cdf_{algorithm_name_metrop.replace('-', '_')}_"
                        f"{external_global_moves_string}_{no_of_sites_string}_{no_of_observations_metrop}_obs_"
                        f"temp_eq_{temperature:.4f}_job_{job_index}.npy", cdfs_of_magnetisation_phase_metrop[job_index])
                if job_index < no_of_jobs_ecmc:
                    np.save(f"{output_directory}/magnetisation_phase_cdf_{algorithm_name_ecmc.replace('-', '_')}_"
                            f"{external_global_moves_string}_{no_of_sites_string}_{no_of_observations_ecmc}_obs_"
                            f"temp_eq_{temperature:.4f}_job_{job_index}.npy",
                            cdfs_of_magnetisation_phase_ecmc[job_index])
                    if subplot_index == 1:
                        inset_axis.plot(*cdfs_of_magnetisation_phase_ecmc[job_index], color=colors[temperature_index],
                                        linestyle=linestyles[job_index], linewidth=2)


def get_cdf_of_magnetisation_phase(sample_directory, temperature, temperature_index, no_of_sites,
                                   no_of_equilibration_sweeps):
    r"""
    Returns (for the sample of the magnetisation phase) the empirical CDF F_n(x), where n is the sample size and F(x)
    is the CDF of the continuous uniform distribution Uniform(-pi, pi).

    The magnetisation phase phi(x; temperature, no_of_sites) is defined via
    m(x; temperature, no_of_sites) = (|| m(x; temperature, no_of_sites) ||, phi(x; temperature, no_of_sites))^t in
    radial coordinates, where m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the
    Cartesian magnetisation, with x_i the position of particle i at the time of observation.

    Parameters
    ----------
    sample_directory : str
        The location of the directory containing the sample and Metropolis acceptance rate.
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.

    Returns
    -------
    numpy.ndarray
        The empirical CDF.  A two-dimensional numpy array (of floats) of size (2, n) representing the empirical CDF.
    """
    return markov_chain_diagnostics.get_cumulative_distribution(sample_getter.get_magnetisation_phase(
        sample_directory, temperature, temperature_index, no_of_sites)[no_of_equilibration_sweeps:])


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: no positional arguments allowed.")
    else:
        main()
