import importlib
import matplotlib
import matplotlib.pyplot as plt
import os
import sys


# Add the directory that contains config_file, sample_getter and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
config_data_getter = importlib.import_module("config_data_getter")
sample_getter = importlib.import_module("sample_getter")
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")


def main(config_file_1, config_file_2):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory_1, output_directory_2, no_of_sites, no_of_equilibration_sweeps_1,
     no_of_equilibration_sweeps_2, temperature) = get_required_config_data_and_error_check(config_file_1, config_file_2)
    beta = 1.0 / temperature
    temperature_directory = f"temp_eq_{temperature:.2f}"

    if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
        sample_1 = sample_getter.get_potential(output_directory_1, temperature_directory, beta, no_of_sites)[
                   no_of_equilibration_sweeps_1:]
        sample_2 = sample_getter.get_potential(output_directory_2, temperature_directory, beta, no_of_sites)[
                   no_of_equilibration_sweeps_2:]
    elif (algorithm_name == "hxy-ecmc" or algorithm_name == "hxy-metropolis" or
          algorithm_name == "hxy-gaussian-noise-metropolis" or algorithm_name == "xy-ecmc" or
          algorithm_name == "xy-metropolis" or algorithm_name == "xy-gaussian-noise-metropolis"):
        sample_1 = sample_getter.get_magnetisation_norm(output_directory_1, temperature_directory, beta, no_of_sites)[
                 no_of_equilibration_sweeps_1:]
        sample_2 = sample_getter.get_magnetisation_norm(output_directory_2, temperature_directory, beta,
                                                        no_of_sites)[no_of_equilibration_sweeps_2:]

    effective_sample_size_1 = markov_chain_diagnostics.get_effective_sample_size(sample_1)
    print(f"Effective sample size (first config file) = {effective_sample_size_1} (from a total sample size of "
          f"{len(sample_1)}).")
    effective_sample_size_2 = markov_chain_diagnostics.get_effective_sample_size(sample_2)
    print(f"Effective sample size (second config file) = {effective_sample_size_2} (from a total sample size of "
          f"{len(sample_2)}).")

    sample_1_cdf = markov_chain_diagnostics.get_cumulative_distribution(sample_1)
    sample_2_cdf = markov_chain_diagnostics.get_cumulative_distribution(sample_2)
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
    (algorithm_name_1, output_directory_1, no_of_sites_1, no_of_equilibration_sweeps_1, initial_temperature_1,
     final_temperature_1, no_of_temperature_increments_1, no_of_jobs_1) = config_data_getter.get_basic_data(
        config_file_1)
    temperature_1 = initial_temperature_1
    (algorithm_name_2, output_directory_2, no_of_sites_2, no_of_equilibration_sweeps_2, initial_temperature_2,
     final_temperature_2, no_of_temperature_increments_2, no_of_jobs_2) = config_data_getter.get_basic_data(
        config_file_2)
    temperature_2 = initial_temperature_2
    if no_of_temperature_increments_1 != 0:
        print("ConfigurationError: in the first configuration file, the value of the no_of_temperature_increments "
              "does not equal 0. In order to compare to single-temperature simulations, this is required.")
        raise SystemExit
    if no_of_temperature_increments_2 != 0:
        print("ConfigurationError: in the second configuration file, the value of the no_of_temperature_increments "
              "does not equal 0. In order to compare to single-temperature simulations, this is required.")
        raise SystemExit
    if no_of_jobs_1 != 1:
        print("ConfigurationError: In the first configuration file, the value of no_of_jobs is not equal to one. Give "
              "configuration files whose value of no_of_jobs is equal to one.")
        raise SystemExit
    if no_of_jobs_2 != 1:
        print("ConfigurationError: In the second configuration file, the value of no_of_jobs is not equal to one. Give"
              " configuration files whose value of no_of_jobs is equal to one.")
        raise SystemExit
    if (((algorithm_name_1 == "hxy-ecmc" or algorithm_name_1 == "hxy-metropolis" or
          algorithm_name_1 == "hxy-gaussian-noise-metropolis") and not
         (algorithm_name_2 == "hxy-ecmc" or algorithm_name_2 == "hxy-metropolis" or
          algorithm_name_2 == "hxy-gaussian-noise-metropolis")) or
            ((algorithm_name_1 == "xy-ecmc" or algorithm_name_1 == "xy-metropolis" or
              algorithm_name_1 == "hxy-gaussian-noise-metropolis") and not
             (algorithm_name_2 == "xy-ecmc" or algorithm_name_2 == "xy-metropolis" or
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
