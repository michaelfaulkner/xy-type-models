from cramer_von_mises import get_cvm_mag_phase_vs_temperature
from setup_scripts import setup_pool
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
    linear_system_sizes = [2 ** (index + 3) for index in range(4)]
    config_files_metrop = [f"config_files/cvm_figure/{value}x{value}_metrop.txt" for value in linear_system_sizes]
    config_files_metrop_local = [f"config_files/cvm_ratio_vs_temp/{value}x{value}_metrop_local.txt" for value in
                                 linear_system_sizes]

    (algorithm_name_metrop, sample_directory_8x8_metrop_all, _, _, no_of_equilibration_sweeps_metrop,
     no_of_observations_metrop, temperatures, _, external_global_moves_string_all, no_of_jobs_metrop, _,
     max_no_of_cpus) = run_script.get_config_data(config_files_metrop[0])
    (_, sample_directory_8x8_metrop_local, _, _, _, _, _, _, external_global_moves_string_local, _, _, _
     ) = run_script.get_config_data(config_files_metrop_local[0])
    main_cvm_output_directory = sample_directory_8x8_metrop_all.replace("/8x8_metrop", "")
    cvm_ratio_output_directory = sample_directory_8x8_metrop_local.replace("/8x8_metrop_local", "")
    pool = setup_pool(no_of_jobs_metrop, max_no_of_cpus)

    figure, axis = plt.subplots(1)
    figure.tight_layout()
    axis.tick_params(which='both', direction='in', length=4, width=2, labelsize=12)
    axis.set_xlabel(r"$1 / (\beta J)$", fontsize=20, labelpad=3)
    axis.set_ylabel(r"$\omega_{\phi_m,n}^{2,{\rm all}} / \omega_{\phi_m,n}^{2,{\rm local}}$", fontsize=20, labelpad=3)
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"]]
    axis.tick_params(which='both', direction='in', width=3)
    axis.tick_params(which='major', length=7, labelsize=18, pad=5)
    axis.tick_params(which='minor', length=4)

    colors = ["black", "red", "blue", "green", "yellow", "cyan"]

    start_time = time.time()
    for system_size_index, length in enumerate(linear_system_sizes):
        try:
            with open(f"{cvm_ratio_output_directory}/magnetisation_phase_cvm_ratio_"
                      f"{algorithm_name_metrop.replace('-', '_')}_{length}x{length}_sites.tsv", "r") as output_file:
                output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                    if not line.startswith('#')]).transpose()
                cvm_ratios = output_file_sans_header[1]
        except IOError:
            cvm_ratio_file = open(f"{cvm_ratio_output_directory}/magnetisation_phase_cvm_ratio_"
                                  f"{algorithm_name_metrop.replace('-', '_')}_{length}x{length}_sites.tsv", "w")
            cvm_ratio_file.write("# temperature".ljust(30) + "CvM ratio".ljust(30) + "CvM ratio error".ljust(30) +
                                 "CvM (all moves)".ljust(30) + "CvM error (all moves)".ljust(30) +
                                 "CvM (local only)".ljust(30) + "CvM error (local only)".ljust(30) + "\n")
            cvm_metrop_alls, cvm_metrop_all_errors = get_cvm_mag_phase_vs_temperature(
                algorithm_name_metrop, main_cvm_output_directory,
                f"{main_cvm_output_directory}/{length}x{length}_metrop", no_of_equilibration_sweeps_metrop,
                no_of_observations_metrop, temperatures, external_global_moves_string_all, no_of_jobs_metrop, pool,
                length)
            cvm_metrop_locals, cvm_metrop_local_errors = get_cvm_mag_phase_vs_temperature(
                algorithm_name_metrop, cvm_ratio_output_directory,
                f"{cvm_ratio_output_directory}/{length}x{length}_metrop_local_moves", no_of_equilibration_sweeps_metrop,
                no_of_observations_metrop, temperatures, external_global_moves_string_local, no_of_jobs_metrop, pool,
                length)
            cvm_ratios = [cvm_metrop_alls[temperature_index] / cvm_metrop_locals[temperature_index]
                          for temperature_index in range(len(temperatures))]
            cvm_ratio_errors = [math.sqrt(
                (cvm_metrop_all_errors[temperature_index] / cvm_metrop_locals[temperature_index]) ** 2 +
                (cvm_metrop_alls[temperature_index] * cvm_metrop_local_errors[temperature_index] /
                 cvm_metrop_locals[temperature_index] ** 2) ** 2) for temperature_index in range(len(temperatures))]
            for temperature_index, temperature in enumerate(temperatures):
                cvm_ratio_file.write(f"{temperature:.14e}".ljust(30) +
                                     f"{cvm_ratios[temperature_index]:.14e}".ljust(30) +
                                     f"{cvm_ratio_errors[temperature_index]:.14e}".ljust(30) +
                                     f"{cvm_metrop_alls[temperature_index]:.14e}".ljust(30) +
                                     f"{cvm_metrop_all_errors[temperature_index]:.14e}".ljust(30) +
                                     f"{cvm_metrop_locals[temperature_index]:.14e}".ljust(30) +
                                     f"{cvm_metrop_local_errors[temperature_index]:.14e}".ljust(30) + "\n")
            cvm_ratio_file.close()
        axis.plot(temperatures, cvm_ratios, marker=".", markersize=5, color=colors[system_size_index],
                  linestyle='dashed', label=fr"$N$ = {length}x{length}")

    legend = axis.legend(loc="upper right", fontsize=10)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(3)
    figure.savefig(f"{cvm_ratio_output_directory}/magnetisation_phase_cvm_ratio_xy_gaussian_noise_metropolis.pdf",
                   bbox_inches="tight")

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()


if __name__ == "__main__":
    if len(sys.argv) != 1:
        raise Exception("InterfaceError: no positional arguments allowed.")
    else:
        main()
