from markov_chain_diagnostics import get_cumulative_distribution, get_effective_sample_size
from sample_getter import get_magnetisation_norm, get_potential
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


def main(config_file_1, config_file_2):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory_1, output_directory_2, no_of_sites, no_of_equilibration_sweeps_1,
     no_of_equilibration_sweeps_2, temperature) = get_required_config_data_and_error_check(config_file_1, config_file_2)

    if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
        sample_1 = get_potential(output_directory_1, temperature, 0, no_of_sites, no_of_equilibration_sweeps_1)
        sample_2 = get_potential(output_directory_2, temperature, 0, no_of_sites, no_of_equilibration_sweeps_2)
    elif (algorithm_name == "hxy-ecmc" or algorithm_name == "hxy-uniform-noise-metropolis" or
          algorithm_name == "hxy-gaussian-noise-metropolis" or algorithm_name == "xy-ecmc" or
          algorithm_name == "xy-uniform-noise-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis"):
        sample_1 = get_magnetisation_norm(output_directory_1, temperature, 0, no_of_sites, no_of_equilibration_sweeps_1)
        sample_2 = get_magnetisation_norm(output_directory_2, temperature, 0, no_of_sites, no_of_equilibration_sweeps_2)
    else:
        sample_1, sample_2 = None, None

    effective_sample_size_1 = get_effective_sample_size(sample_1)
    print(f"Effective sample size (first config file) = {effective_sample_size_1} (from a total sample size of "
          f"{len(sample_1)}).")
    effective_sample_size_2 = get_effective_sample_size(sample_2)
    print(f"Effective sample size (second config file) = {effective_sample_size_2} (from a total sample size of "
          f"{len(sample_2)}).")

    sample_1_cdf = get_cumulative_distribution(sample_1)
    sample_2_cdf = get_cumulative_distribution(sample_2)
    plt.plot(sample_1_cdf[0], sample_1_cdf[1], color="r", linewidth=4, linestyle="-", label="first config file")
    plt.plot(sample_2_cdf[0], sample_2_cdf[1], color="k", linewidth=2, linestyle="-", label="second config file")

    plt.xlabel(r"$x$", fontsize=15, labelpad=10)
    plt.ylabel(r"$ \mathbb{P} \left( X < x \right)$", fontsize=15, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
    legend = plt.legend(loc="upper left", fontsize=10)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(1.5)
    plt.tight_layout()
    plt.show()


def get_required_config_data_and_error_check(config_file_1, config_file_2):
    (algorithm_name_1, output_directory_1, no_of_sites_1, _, no_of_equilibration_sweeps_1, _, temperatures_1, _, _,
     no_of_runs_1, initial_run_index_1, max_no_of_cpus_1) = run_script.get_config_data(config_file_1)
    temperature_1 = temperatures_1[0]
    (algorithm_name_2, output_directory_2, no_of_sites_2, _, no_of_equilibration_sweeps_2, _, temperatures_2, _, _,
     no_of_runs_2, initial_run_index_2, max_no_of_cpus_2) = run_script.get_config_data(config_file_2)
    temperature_2 = temperatures_2[0]
    if len(temperatures_1) != 1:
        print("ConfigurationError: In order to compare single-temperature simulations, set the value of "
              "no_of_temperature_increments equal to 0 in both configuration files.")
        raise SystemExit
    if len(temperatures_2) != 1:
        print("ConfigurationError: In order to compare single-temperature simulations, set the value of "
              "no_of_temperature_increments equal to 0 in both configuration files.")
        raise SystemExit
    if no_of_runs_1 != 1 or initial_run_index_1 != 0:
        print("ConfigurationError: Set the value of no_of_runs to 1 and the value of initial_run_index to 0 in both "
              "configuration files.")
        raise SystemExit
    if no_of_runs_2 != 1 or initial_run_index_2 != 0:
        print("ConfigurationError: Set the value of no_of_runs to 1 and the value of initial_run_index to 0 in both "
              "configuration files.")
        raise SystemExit
    if (((algorithm_name_1 == "hxy-ecmc" or algorithm_name_1 == "hxy-uniform-noise-metropolis" or
          algorithm_name_1 == "hxy-gaussian-noise-metropolis") and not
         (algorithm_name_2 == "hxy-ecmc" or algorithm_name_2 == "hxy-uniform-noise-metropolis" or
          algorithm_name_2 == "hxy-gaussian-noise-metropolis")) or
            ((algorithm_name_1 == "xy-ecmc" or algorithm_name_1 == "xy-uniform-noise-metropolis" or
              algorithm_name_1 == "hxy-gaussian-noise-metropolis") and not
             (algorithm_name_2 == "xy-ecmc" or algorithm_name_2 == "xy-uniform-noise-metropolis" or
              algorithm_name_2 == "xy-gaussian-noise-metropolis")) or
            (algorithm_name_1 == "elementary-electrolyte" and algorithm_name_2 != "elementary-electrolyte") or
            (algorithm_name_1 == "multivalued-electrolyte" and algorithm_name_2 != "multivalued-electrolyte")):
        print("ConfigurationError: give the same model in each configuration file.")
        raise SystemExit
    if temperature_1 != temperature_2:
        print("ConfigurationError: give the same initial_temperature in each configuration file.")
        raise SystemExit
    if no_of_sites_1 != no_of_sites_2:
        print("ConfigurationError: give the same integer_lattice_length in each configuration file.")
        raise SystemExit
    return (algorithm_name_1, output_directory_1, output_directory_2, no_of_sites_1, no_of_equilibration_sweeps_1,
            no_of_equilibration_sweeps_2, temperature_1)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise Exception("InterfaceError: Two positional arguments required - give the locations of the two "
                        "configuration files.")
    main(sys.argv[1], sys.argv[2])
