from sample_getter import get_cartesian_magnetisation
from setup_scripts import reverse_enumerate
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


def main(symmetry_breaking_paper=True):
    """For 'General symmetry breaking at a topological phase transition,' set symmetry_breaking_paper=True;
        for 'Sampling algorithms in statistical physics...', set symmetry_breaking_paper=False."""
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    if symmetry_breaking_paper:
        config_file_16x16_metrop = "config_files/mag_evolution_figs/16x16_metrop.txt"
        config_file_64x64_metrop = "config_files/mag_evolution_figs/64x64_metrop.txt"
        config_file_256x256_metrop = "config_files/mag_evolution_figs/256x256_metrop.txt"
        config_file_16x16_ecmc = "config_files/mag_evolution_figs/16x16_ecmc.txt"
        config_file_64x64_ecmc = "config_files/mag_evolution_figs/64x64_ecmc.txt"
        config_file_256x256_ecmc = "config_files/mag_evolution_figs/256x256_ecmc.txt"
        config_file_256x256_metrop_with_twists = "config_files/mag_evolution_figs/256x256_metrop_with_twists.txt"
    else:
        config_file_16x16_metrop = "config_files/sampling_algos_xy_figs/16x16_metrop.txt"
        config_file_64x64_metrop = "config_files/sampling_algos_xy_figs/64x64_metrop.txt"
        config_file_256x256_metrop = "config_files/sampling_algos_xy_figs/256x256_metrop.txt"
        config_file_16x16_ecmc = "config_files/sampling_algos_xy_figs/16x16_ecmc.txt"
        config_file_64x64_ecmc = "config_files/sampling_algos_xy_figs/64x64_ecmc.txt"
        config_file_256x256_ecmc = "config_files/sampling_algos_xy_figs/256x256_ecmc.txt"
        config_file_256x256_metrop_with_twists = None

    (algorithm_name_metrop, sample_directory_16x16_metrop, no_of_sites_16x16, no_of_sites_string_16x16,
     no_of_equilibration_sweeps_metrop, no_of_samples_metrop, temperatures, use_external_global_moves_16x16_metrop,
     external_global_moves_string_16x16_metrop, no_of_runs_metrop, _, max_no_of_cpus) = run_script.get_config_data(
        config_file_16x16_metrop)
    (_, sample_directory_64x64_metrop, no_of_sites_64x64, no_of_sites_string_64x64, _, _, _,
     use_external_global_moves_64x64_metrop, external_global_moves_string_64x64_metrop, _, _, _
     ) = run_script.get_config_data(config_file_64x64_metrop)
    (_, sample_directory_256x256_metrop, no_of_sites_256x256, no_of_sites_string_256x256, _, _, _,
     use_external_global_moves_256x256_metrop, external_global_moves_string_256x256_metrop, _, _, _
     ) = run_script.get_config_data(config_file_256x256_metrop)
    (algorithm_name_ecmc, sample_directory_16x16_ecmc, _, _, no_of_equilibration_sweeps_ecmc, no_of_samples_ecmc,
     _, use_external_global_moves_ecmc, external_global_moves_string_ecmc, no_of_runs_ecmc, _, _
     ) = run_script.get_config_data(config_file_16x16_ecmc)
    sample_directory_64x64_ecmc = run_script.get_config_data(config_file_64x64_ecmc)[1]
    sample_directory_256x256_ecmc = run_script.get_config_data(config_file_256x256_ecmc)[1]

    output_directory = sample_directory_16x16_metrop.replace("/16x16_metrop", "")
    alphabetic_label_16x16_metrop = "a"
    alphabetic_label_64x64_metrop = "b"
    alphabetic_label_256x256_metrop = "c"
    alphabetic_label_16x16_ecmc = "d"
    alphabetic_label_64x64_ecmc = "e"
    alphabetic_label_256x256_ecmc = "f"

    start_time = time.time()
    # for run_index in range(no_of_runs_metrop_with_twists):
    figure, axes = plt.subplots(2, 3, figsize=(30, 20))
    make_subplot(axes[0, 0], algorithm_name_metrop, output_directory, sample_directory_16x16_metrop,
                 no_of_sites_16x16, no_of_sites_string_16x16, no_of_equilibration_sweeps_metrop,
                 no_of_samples_metrop, temperatures, use_external_global_moves_16x16_metrop,
                 external_global_moves_string_16x16_metrop, alphabetic_label_16x16_metrop, None)
    make_subplot(axes[0, 1], algorithm_name_metrop, output_directory, sample_directory_64x64_metrop,
                 no_of_sites_64x64, no_of_sites_string_64x64, no_of_equilibration_sweeps_metrop,
                 no_of_samples_metrop, temperatures, use_external_global_moves_64x64_metrop,
                 external_global_moves_string_64x64_metrop, alphabetic_label_64x64_metrop, None)
    make_subplot(axes[0, 2], algorithm_name_metrop, output_directory, sample_directory_256x256_metrop,
                 no_of_sites_256x256, no_of_sites_string_256x256, no_of_equilibration_sweeps_metrop,
                 no_of_samples_metrop, temperatures, use_external_global_moves_256x256_metrop,
                 external_global_moves_string_256x256_metrop, alphabetic_label_256x256_metrop, None)
    make_subplot(axes[1, 0], algorithm_name_ecmc, output_directory, sample_directory_16x16_ecmc,
                 no_of_sites_16x16, no_of_sites_string_16x16, no_of_equilibration_sweeps_ecmc,
                 no_of_samples_ecmc, temperatures, use_external_global_moves_ecmc,
                 external_global_moves_string_ecmc, alphabetic_label_16x16_ecmc, None)
    make_subplot(axes[1, 1], algorithm_name_ecmc, output_directory, sample_directory_64x64_ecmc,
                 no_of_sites_64x64, no_of_sites_string_64x64, no_of_equilibration_sweeps_ecmc,
                 no_of_samples_ecmc, temperatures, use_external_global_moves_ecmc,
                 external_global_moves_string_ecmc, alphabetic_label_64x64_ecmc, None)
    make_subplot(axes[1, 2], algorithm_name_ecmc, output_directory, sample_directory_256x256_ecmc,
                 no_of_sites_256x256, no_of_sites_string_256x256, no_of_equilibration_sweeps_ecmc,
                 no_of_samples_ecmc, temperatures, use_external_global_moves_ecmc,
                 external_global_moves_string_ecmc, alphabetic_label_256x256_ecmc, None)

    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
    plt.subplots_adjust(wspace=0.125, hspace=0.05)
    figure.add_artist(
        plt.Line2D([0.115, 0.91], [0.495, 0.495], transform=figure.transFigure, color="black", linewidth=4.5))
    figure.add_artist(
        plt.Line2D([0.378, 0.378], [0.11, 0.88], transform=figure.transFigure, color="black", linewidth=4.5))
    figure.add_artist(
        plt.Line2D([0.6465, 0.6465], [0.11, 0.88], transform=figure.transFigure, color="black", linewidth=4.5))

    figure.savefig(f"{output_directory}/magnetisation_evolution_xy_model_metropolis_with_gaussian_noise_and_ecmc_"
                   f"{no_of_samples_metrop}_metrop_obs_{no_of_samples_ecmc}_ecmc_obs.png",
                   bbox_inches="tight")

    if symmetry_breaking_paper:
        (_, sample_directory_256x256_metrop_with_twists, _, _, _, no_of_samples_256x256_metrop_with_twists,
         temperatures_256x256_metrop_with_twists, use_external_global_moves_256x256_metrop_with_twists,
         external_global_moves_string_256x256_metrop_with_twists, no_of_runs_metrop_with_twists, _, _
         ) = run_script.get_config_data(config_file_256x256_metrop_with_twists)
        for run_index in range(no_of_runs_metrop_with_twists):
            figure.clear()
            figure, axis = plt.subplots(1)
            make_subplot(axis, algorithm_name_metrop, output_directory, sample_directory_256x256_metrop_with_twists,
                         no_of_sites_256x256, no_of_sites_string_256x256, no_of_equilibration_sweeps_metrop,
                         no_of_samples_256x256_metrop_with_twists, temperatures_256x256_metrop_with_twists,
                         use_external_global_moves_256x256_metrop_with_twists,
                         external_global_moves_string_256x256_metrop_with_twists, None, run_index)
            figure.savefig(
                f"{output_directory}/magnetisation_evolution_xy_model_metropolis_with_gaussian_noise_and_global_twists_"
                f"{no_of_samples_256x256_metrop_with_twists}_metrop_obs_run_{run_index}.png", bbox_inches="tight")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


def make_subplot(axis, algorithm_name, output_directory, sample_directory, no_of_sites, no_of_sites_string,
                 no_of_equilibration_sweeps, no_of_samples, temperatures, use_external_global_moves,
                 external_global_moves_string, alphabetic_label, run_index):
    axis.axis('square')
    if alphabetic_label is None:
        axis.text(-0.05, 0.96, "(a)", fontsize=22, transform=axis.transAxes)
        """left old code in comment below as the bbox=dict() part may be useful"""
        # axis.text(0.76125, 0.78, fr"$N = {int(no_of_sites ** 0.5)} \! \times \! {int(no_of_sites ** 0.5)}$",
        #           fontsize=14, transform=axis.transAxes,
        #           bbox=dict(facecolor='none', edgecolor='black', linewidth=3, boxstyle='round, pad=0.5'))
        axis.text(0.69, 0.94, fr"$N = {int(no_of_sites ** 0.5)} \! \times \! {int(no_of_sites ** 0.5)}$",
                  fontsize=15, transform=axis.transAxes,
                  bbox=dict(facecolor='none', edgecolor='black', linewidth=3.0, boxstyle='round, pad=0.5'))
        """add text box to pad the bottom of plot - aids formatting when placed next to twist-prob fig"""
        axis.text(1.0, -0.12, " ", fontsize=16, transform=axis.transAxes,
                  bbox=dict(facecolor='none', edgecolor='none', linewidth=1.0, boxstyle='round, pad=0.5'))
        axis.tick_params(which='both', length=5.333, width=3)
        axis.tick_params(which='minor', length=4)
    else:
        axis.text(-0.046, 0.9675, f"({alphabetic_label})", fontsize=40, transform=axis.transAxes)
        if alphabetic_label == "c" or alphabetic_label == "f":
            axis.text(0.692, 0.965, fr"$N = {int(no_of_sites ** 0.5)} \! \times \! {int(no_of_sites ** 0.5)}$",
                      fontsize=25, transform=axis.transAxes,
                      bbox=dict(facecolor='none', edgecolor='black', linewidth=4.5, boxstyle='round, pad=0.5'))
        else:
            axis.text(0.75, 0.965, fr"$N = {int(no_of_sites ** 0.5)} \! \times \! {int(no_of_sites ** 0.5)}$",
                      fontsize=25, transform=axis.transAxes,
                      bbox=dict(facecolor='none', edgecolor='black', linewidth=4.5, boxstyle='round, pad=0.5'))
        axis.tick_params(which='both', length=8, width=4.5)
        axis.tick_params(which='minor', length=6)
    if alphabetic_label is None:
        axis.tick_params(which='major', labelsize=21, pad=2.5)
    elif alphabetic_label == "a":
        axis.tick_params(which='major', labelsize=24, pad=10)
    else:
        axis.tick_params(which='major', labelsize=30, pad=10)

    # create four-quadrant axes
    if alphabetic_label is None:
        axis.spines['bottom'].set_position('zero')
        axis.spines['left'].set_linewidth(3), axis.spines['bottom'].set_linewidth(3)
    else:
        axis.spines['bottom'].set_position('center')
        axis.spines['left'].set_linewidth(4.5), axis.spines['bottom'].set_linewidth(4.5)
    axis.spines['left'].set_position('center')
    axis.spines['top'].set_color('none'), axis.spines['right'].set_color('none')
    axis.xaxis.set_ticks_position('bottom'), axis.yaxis.set_ticks_position('left')

    minor_ticks = np.arange(-0.75, 1.0, 0.25)
    axis.xaxis.set_major_formatter('{x:.1f}')
    axis.set_xlim([-1.0, 1.0])
    axis.set_xticks([-1.0, 1.0])
    axis.set_xticks(minor_ticks, minor=True)
    axis.yaxis.set_major_formatter('{x:.1f}')
    if alphabetic_label is None:
        axis.set_ylim([-0.72, 1.0])
        axis.set_yticks(np.arange(-0.5, 1.0, 0.25), minor=True)
        axis.set_yticks([1.0])
        axis.set_xlabel(r"$m_x$", fontsize=21), axis.set_ylabel(r"$m_y$", fontsize=21, rotation="horizontal")
        axis.xaxis.set_label_coords(1.005, 0.505), axis.yaxis.set_label_coords(0.57, 0.95)
    else:
        axis.set_ylim([-1.0, 1.0])
        axis.set_yticks(minor_ticks, minor=True)
        axis.set_yticks([-1.0, 1.0])
        axis.set_xlabel(r"$m_x$", fontsize=35), axis.set_ylabel(r"$m_y$", fontsize=35, rotation="horizontal")
        axis.xaxis.set_label_coords(1.005, 0.575), axis.yaxis.set_label_coords(0.565, 0.96)

    colors = ["black", "red"]
    for temperature_index, temperature in reverse_enumerate(temperatures):
        print(f"Temperature = {temperature:.4f}")
        try:
            if use_external_global_moves:
                cartesian_magnetisation = np.load(
                    f"{output_directory}/cartesian_magnetisation_sample_{algorithm_name.replace('-', '_')}_"
                    f"{external_global_moves_string}_{no_of_sites_string}_{no_of_samples}_obs_temp_eq_"
                    f"{temperature:.4f}_run_{run_index}.npy")
            else:
                cartesian_magnetisation = np.load(
                    f"{output_directory}/cartesian_magnetisation_sample_{algorithm_name.replace('-', '_')}_"
                    f"{external_global_moves_string}_{no_of_sites_string}_{no_of_samples}_obs_temp_eq_"
                    f"{temperature:.4f}.npy")
        except IOError:
            if use_external_global_moves:
                cartesian_magnetisation = get_cartesian_magnetisation(f"{sample_directory}/run_{run_index}",
                                                                      temperature, temperature_index, no_of_sites,
                                                                      no_of_equilibration_sweeps)
                np.save(f"{output_directory}/cartesian_magnetisation_sample_{algorithm_name.replace('-', '_')}_"
                        f"{external_global_moves_string}_{no_of_sites_string}_{no_of_samples}_obs_temp_eq_"
                        f"{temperature:.4f}_run_{run_index}.npy", cartesian_magnetisation)
            else:
                cartesian_magnetisation = get_cartesian_magnetisation(sample_directory, temperature, temperature_index,
                                                                      no_of_sites, no_of_equilibration_sweeps)
                np.save(f"{output_directory}/cartesian_magnetisation_sample_{algorithm_name.replace('-', '_')}_"
                        f"{external_global_moves_string}_{no_of_sites_string}_{no_of_samples}_obs_temp_eq_"
                        f"{temperature:.4f}.npy", cartesian_magnetisation)
        if alphabetic_label is None:
            rotation = 0.0  # previously set rotation = -0.15
            axis.plot(cartesian_magnetisation[:, 0] * math.cos(rotation) +
                      cartesian_magnetisation[:, 1] * math.sin(rotation),
                      cartesian_magnetisation[:, 1] * math.cos(rotation) -
                      cartesian_magnetisation[:, 0] * math.sin(rotation), linestyle="solid", linewidth=1.0,
                      color=colors[temperature_index], label=fr"$1 / (\beta J)$ = {temperature:.2f}")
        else:
            axis.plot(cartesian_magnetisation[:, 0], cartesian_magnetisation[:, 1], linestyle="solid", linewidth=1.0,
                      color=colors[temperature_index], label=fr"$1 / (\beta J)$ = {temperature:.2f}")
    handles, labels = axis.get_legend_handles_labels()
    if alphabetic_label is None:
        legend = axis.legend(reversed(handles), reversed(labels), loc='lower right', bbox_to_anchor=(1.07, -0.0325),
                             fontsize=15)
    elif alphabetic_label == "a":
        legend = axis.legend(reversed(handles), reversed(labels), loc='lower right', bbox_to_anchor=(1.06, -0.05),
                             fontsize=17)
    elif alphabetic_label == "d":
        legend = axis.legend(reversed(handles), reversed(labels), loc='lower right', bbox_to_anchor=(1.06, -0.05),
                             fontsize=18)
    else:
        legend = axis.legend(reversed(handles), reversed(labels), loc='lower right', bbox_to_anchor=(1.06, -0.05),
                             fontsize=20)
    legend.get_frame().set_edgecolor("black")
    if alphabetic_label is None:
        legend.get_frame().set_lw(3)
    else:
        legend.get_frame().set_lw(4.5)
    [line.set_linewidth(2) for line in legend.get_lines()]


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("InterfaceError: At most one positional argument permitted.  None are required but you may "
                        "provide symmetry_breaking_paper, a Boolean (default value is True, corresponding to the "
                        "symmetry-breaking paper (see docstring of main()).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be symmetry_breaking_paper, a Boolean (default value is "
              "True, corresponding to the symmetry-breaking paper (see docstring of main()).")
        if sys.argv[1] == "True":
            main(True)
        elif sys.argv[1] == "False":
            main(False)
        else:
            Exception("InterfaceError: The provided value of symmetry_breaking_paper is neither True nor False (see "
                      "docstring of main()).")
    else:
        print("No positional arguments provided.  None are required but you may provide symmetry_breaking_paper, a "
              "Boolean (default value is True, corresponding to the symmetry-breaking paper (see docstring of main()).")
        main()
