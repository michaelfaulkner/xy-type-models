import importlib
import matplotlib
import matplotlib.pyplot as plt
import os
import sys

# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
import numpy as np

this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
setup_scripts = importlib.import_module("setup_scripts")
sample_getter = importlib.import_module("sample_getter")


def main(config_file, observable_string, length_of_trace_plot=1000, number_of_histogram_bins=1000):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = setup_scripts.get_config_data(config_file)
    setup_scripts.check_for_observable_error(algorithm_name, observable_string)
    (temperature, magnitude_of_temperature_increments) = setup_scripts.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)
    if no_of_jobs != 1:
        sample_directory = output_directory + f"/job_1"
    else:
        sample_directory = output_directory

    for i in range(no_of_temperature_increments + 1):
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"
        get_sample_method = getattr(sample_getter, "get_" + observable_string)
        sample = np.atleast_2d(get_sample_method(sample_directory, temperature_directory, beta, no_of_sites))

        plt.xlabel(r"time, $t$", fontsize=15, labelpad=10)
        plt.ylabel(r"$ x \left( t \right)$", fontsize=15, labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
        plt.plot(sample[0, no_of_equilibration_sweeps:no_of_equilibration_sweeps + length_of_trace_plot], color="k",
                 linewidth=1, linestyle="-")
        plt.tight_layout()
        plt.savefig(f"{output_directory}/{observable_string}_vs_time_temp_eq_{temperature:.2f}_"
                    f"{int(no_of_sites ** 0.5)}_{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf",
                    bbox_inches="tight")
        plt.clf()

        plt.xlabel(r"$x$", fontsize=15, labelpad=10)
        plt.ylabel(r"$ \mathbb{P} \left( |x| < |X| < |x| + |dx| \right)$", fontsize=15, labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
        plt.hist(sample[0, no_of_equilibration_sweeps:] - np.mean(sample[0, no_of_equilibration_sweeps:]),
                 bins=number_of_histogram_bins)
        plt.savefig(f"{output_directory}/{observable_string}_histogram_temp_eq_{temperature:.2f}_"
                    f"{int(no_of_sites ** 0.5)}_{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf",
                    bbox_inches="tight")

        temperature = temperature - magnitude_of_temperature_increments


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose power trispectrum you wish to estimate.  In addition, you "
                        "may provide length_of_trace_plot and number_of_histogram_bins (both default values are 1000) "
                        "in the third and fourth positions (respectively).")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  In addition, you may provide length_of_trace_plot and "
              "number_of_histogram_bins (both default values are 1000) in the third and fourth positions "
              "(respectively).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The third must be length_of_trace_plot.  In addition, you may "
              "provide number_of_histogram_bins (default value is 1000) in the fourth position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The third / fourth must be length_of_trace_plot / "
              "number_of_histogram_bins.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
