import importlib
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import time

# Add the directory that contains config_file and markov_chain_diagnostics to sys.path
import numpy as np

this_directory = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, this_directory)
setup_scripts = importlib.import_module("setup_scripts")
sample_getter = importlib.import_module("sample_getter")


def main(config_file, observable_string, max_physical_time=100.0, number_of_histogram_bins=1000):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = setup_scripts.get_config_data(config_file)
    setup_scripts.check_for_observable_error(algorithm_name, observable_string)
    (temperature, magnitude_of_temperature_increments) = setup_scripts.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)
    if no_of_jobs != 1:
        sample_directory = f"{output_directory}/job_1"
    else:
        sample_directory = output_directory

    start_time = time.time()
    for i in range(no_of_temperature_increments + 1):
        print(f"Temperature = {temperature:.2f}")
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"
        physical_time_step = 1.0 / sample_getter.get_sampling_frequency(sample_directory, None, temperature_directory)
        length_of_trace_plot = int(max_physical_time / physical_time_step)
        get_sample_method = getattr(sample_getter, "get_" + observable_string)
        sample = np.atleast_2d(get_sample_method(sample_directory, temperature_directory, beta, no_of_sites))
        if len(sample) > 1:
            sample = sample.transpose()

        plt.xlabel(r"time, $t$ ($s$)", fontsize=15, labelpad=10)
        if len(sample) > 1:
            plt.ylabel(r"$ x_1 \left( t \right) - \bar{x_1}$", fontsize=15, labelpad=10)
        else:
            plt.ylabel(r"$ x \left( t \right) - \bar{x}$", fontsize=15, labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
        plt.plot(np.arange(length_of_trace_plot) * physical_time_step,
                 sample[0, no_of_equilibration_sweeps:no_of_equilibration_sweeps + length_of_trace_plot]
                 - np.mean(sample[0, no_of_equilibration_sweeps:no_of_equilibration_sweeps + length_of_trace_plot]),
                 color="k", linewidth=1, linestyle="-")
        plt.tight_layout()
        plt.savefig(f"{output_directory}/{observable_string}_vs_time_temp_eq_{temperature:.2f}_"
                    f"{int(no_of_sites ** 0.5)}_{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf",
                    bbox_inches="tight")
        plt.clf()

        if len(sample) > 1:
            plt.xlabel(r"$x_1 - \bar{x_1}$", fontsize=15, labelpad=10)
            plt.ylabel(r"$\pi \left( x_1 - \bar{x_1} \right)$ / const", fontsize=15, labelpad=10)
        else:
            plt.xlabel(r"$x - \bar{x}$", fontsize=15, labelpad=10)
            plt.ylabel(r"$\pi \left( x - \bar{x} \right)$ / const", fontsize=15, labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
        plt.hist(sample[0, no_of_equilibration_sweeps:] - np.mean(sample[0, no_of_equilibration_sweeps:]),
                 bins=number_of_histogram_bins, density=True)
        plt.savefig(f"{output_directory}/{observable_string}_histogram_temp_eq_{temperature:.2f}_"
                    f"{int(no_of_sites ** 0.5)}_{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf",
                    bbox_inches="tight")
        plt.clf()

        temperature -= magnitude_of_temperature_increments
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose power trispectrum you wish to estimate.  In addition, you "
                        "may provide max_physical_time (default value is 100.0) and number_of_histogram_bins (default "
                        "value is 1000) in the third and fourth positions (respectively).")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  In addition, you may provide max_physical_time (default value is "
              "100.0) and number_of_histogram_bins (default value is 1000) in the third and fourth positions "
              "(respectively).")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The third must be max_physical_time.  In addition, you may "
              "provide number_of_histogram_bins (default value is 1000) in the fourth position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The third / fourth must be max_physical_time / "
              "number_of_histogram_bins.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
