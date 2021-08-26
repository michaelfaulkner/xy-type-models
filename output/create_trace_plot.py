import importlib
import matplotlib
import matplotlib.pyplot as plt
import os
import sys

# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
config_data_getter = importlib.import_module("config_data_getter")
sample_getter = importlib.import_module("sample_getter")


def main(config_file, observable_string):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = config_data_getter.get_basic_data(config_file)
    config_data_getter.check_for_observable_error(algorithm_name, observable_string)
    if no_of_jobs != 1:
        raise Exception("ConfigurationError: Give a configuration file whose value of no_of_jobs is equal to one.")
    temperature = final_temperature

    beta = 1.0 / temperature
    temperature_directory = f"temp_eq_{temperature:.2f}"
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
    sample = get_sample_method(output_directory, temperature_directory, beta, no_of_sites)

    plt.xlabel(r"time, $t$", fontsize=15, labelpad=10)
    plt.ylabel(r"$ x \left( t \right)$", fontsize=15, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
    plt.plot(sample, color="k", linewidth=1, linestyle="-")
    plt.tight_layout()
    plt.show()
    plt.clf()
    plt.plot(sample[0:no_of_equilibration_sweeps], color="k", linewidth=1, linestyle="-")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose trace plot you wish to create.")
    main(sys.argv[1], sys.argv[2])
