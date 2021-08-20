import importlib
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys


# Add the directory that contains config_file, sample_getter and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_modules = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_modules)
config_data_getter = importlib.import_module("config_data_getter")
sample_getter = importlib.import_module("sample_getter")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")


def main(config_file):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = config_data_getter.get_basic_data(config_file)
    if no_of_jobs != 1:
        print("ConfigurationError: Give a configuration file whose value of no_of_jobs is equal to one.")
        raise SystemExit
    temperature = initial_temperature
    beta = 1.0 / temperature
    temperature_directory = f"temp_eq_{temperature:.2f}"

    if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
        sample = sample_getter.get_potential(output_directory, temperature_directory, beta, no_of_sites)[
                 no_of_equilibration_sweeps:]
    elif (algorithm_name == "hxy-ecmc" or algorithm_name == "hxy-metropolis" or
          algorithm_name == "hxy-gaussian-noise-metropolis" or algorithm_name == "xy-ecmc" or
          algorithm_name == "xy-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis"):
        sample = sample_getter.get_magnetisation_norm(output_directory, temperature_directory, beta,
                                                      no_of_sites)[no_of_equilibration_sweeps:]

    if algorithm_name == "elementary-electrolyte":
        reference_sample = np.loadtxt("output/convergence_tests/electrolyte/elementary/elementary_electrolyte_8x8_"
                                      "sites_temp_1_point_5_potential_reference_sample.dat",
                                      dtype=float)
    if algorithm_name == "multivalued-electrolyte":
        reference_sample = np.loadtxt("output/convergence_tests/electrolyte/multivalued/multivalued_electrolyte_8x8_"
                                      "sites_temp_1_point_5_potential_reference_sample.dat",
                                      dtype=float)
    if (algorithm_name == "hxy-ecmc" or algorithm_name == "hxy-metropolis" or
            algorithm_name == "hxy-gaussian-noise-metropolis"):
        reference_sample = np.loadtxt("output/convergence_tests/hxy/hxy_8x8_sites_temp_1_point_3_magnetisation_norm_"
                                      "reference_sample.dat", dtype=float)
    if (algorithm_name == "xy-ecmc" or algorithm_name == "xy-metropolis" or
            algorithm_name == "xy-gaussian-noise-metropolis"):
        reference_sample = np.loadtxt("output/convergence_tests/xy/xy_8x8_sites_temp_0_point_8_magnetisation_norm_"
                                      "reference_sample.dat", dtype=float)

    effective_sample_size = markov_chain_diagnostics.get_effective_sample_size(sample)
    print(f"Effective sample size = {effective_sample_size} (from a total sample size of {len(sample)}).")

    reference_cdf = markov_chain_diagnostics.get_cumulative_distribution(reference_sample)
    sample_cdf = markov_chain_diagnostics.get_cumulative_distribution(sample)

    plt.plot(reference_cdf[0], reference_cdf[1], color="r", linewidth=4, linestyle="-", label="reference data")
    plt.plot(sample_cdf[0], sample_cdf[1], color="k", linewidth=2, linestyle="-", label="xy-type-models data")

    plt.xlabel(r"$x$", fontsize=15, labelpad=10)
    plt.ylabel(r"$ \mathbb{P} \left( X < x \right)$", fontsize=15, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
    legend = plt.legend(loc="upper left", fontsize=10)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(1.5)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise Exception("InterfaceError: One positional argument required - give the location of the configuration "
                        "file.")
    main(sys.argv[1])
