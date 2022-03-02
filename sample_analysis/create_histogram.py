import importlib
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import time

# import additional modules; have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
setup_scripts = importlib.import_module("setup_scripts")
sample_getter = importlib.import_module("sample_getter")
run_script = importlib.import_module("run")


def main(config_file, observable_string, no_of_histogram_bins=100):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, _, initial_temperature,
     final_temperature, no_of_temperature_increments, _, external_global_moves_string, no_of_jobs,
     max_no_of_cpus) = run_script.get_config_data(config_file)
    setup_scripts.check_for_observable_error(algorithm_name, observable_string)
    (temperature, magnitude_of_temperature_increments) = setup_scripts.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)
    temperature = initial_temperature
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
    sample_is_one_dimensional = setup_scripts.get_sample_is_one_dimensional(observable_string)

    start_time = time.time()
    for temp_index in range(no_of_temperature_increments + 1):
        print(f"Temperature = {temperature:.2f}")

        if sample_is_one_dimensional:
            plt.xlabel(r"$x - \bar{x}$", fontsize=15, labelpad=10)
            plt.ylabel(r"$\pi \left( x \right)$", fontsize=15, labelpad=10)
        else:
            plt.xlabel(r"$x_1 - \bar{x_1}$", fontsize=15, labelpad=10)
            plt.ylabel(r"$\pi \left( x_1 \right)$", fontsize=15, labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
        if no_of_jobs != 1:
            sample = get_sample_method(f"{output_directory}/job_1", temperature, temp_index, no_of_sites)
        else:
            sample = get_sample_method(output_directory, temperature, temp_index, no_of_sites)
        plt.hist(sample[no_of_equilibration_sweeps:] - np.mean(sample[no_of_equilibration_sweeps:]),
                 bins=no_of_histogram_bins, density=True)
        plt.savefig(f"{output_directory}/{observable_string}_histogram_{algorithm_name.replace('-', '_')}_"
                    f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites_temp_eq_"
                    f"{temperature:.2f}.pdf", bbox_inches="tight")
        plt.clf()
        temperature += magnitude_of_temperature_increments
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose trace plot and histogram you wish to estimate (in the first"
                        " and second positions, respectively).  In addition, you may provide no_of_histogram_bins "
                        "(default value is 100) in the third position.")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  The first / second must be the location of the configuration file / "
              "the string of the observable whose trace plot and histogram you wish to estimate.  In addition, you may "
              "provide no_of_histogram_bins (default value is 100) in the third position.")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The first / second / third must be the location of the "
              "configuration file / the string of the observable whose trace plot and histogram you wish to estimate / "
              "no_of_histogram_bins.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
