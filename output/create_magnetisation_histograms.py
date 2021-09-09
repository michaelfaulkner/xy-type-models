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


def main(config_file, number_of_histogram_bins=1000):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, no_of_jobs) = setup_scripts.get_config_data(config_file)
    if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
        print("ConfigurationError: The configuration file corresponds to a Maggs-electrolyte model but this script "
              "requires the XY of HXY model.")
        raise SystemExit
    (temperature, magnitude_of_temperature_increments) = setup_scripts.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)

    if no_of_jobs == 1:
        sample_directories = np.atleast_1d([output_directory])
    elif no_of_jobs < 4:
        sample_directories = np.array([f"{output_directory}/job_{job_number + 1}" for job_number in range(no_of_jobs)])
    else:
        sample_directories = np.array([f"{output_directory}/job_{job_number + 1}" for job_number in range(4)])

    start_time = time.time()
    for i in range(no_of_temperature_increments + 1):
        print(f"Temperature = {temperature:.2f}")
        beta = 1.0 / temperature
        temperature_directory = f"temp_eq_{temperature:.2f}"
        for job_index, sample_directory in enumerate(sample_directories):
            print(f"job {job_index + 1}")
            cartesian_magnetisation_sample = sample_getter.get_cartesian_magnetisation(sample_directory,
                                                                                       temperature_directory, beta,
                                                                                       no_of_sites).transpose()
            magnetisation_norm_sample = sample_getter.get_magnetisation_norm(sample_directory, temperature_directory,
                                                                             beta, no_of_sites)
            magnetisation_phase_sample = sample_getter.get_magnetisation_phase(sample_directory, temperature_directory,
                                                                               beta, no_of_sites)

            figure, axis = plt.subplots(3, 2, figsize=(10, 10))
            axis[2, 0].set_xlabel(r"$x$", fontsize=15, labelpad=10)
            axis[2, 1].set_xlabel(r"$x$", fontsize=15, labelpad=10)
            axis[0, 0].set_ylabel(r"$\pi \left( x = m_{x / y} \right)$", fontsize=15, labelpad=10)
            axis[1, 0].set_ylabel(r"$\pi \left( x = |m_{x / y}| \right)$", fontsize=15, labelpad=10)
            axis[2, 0].set_ylabel(r"$\pi \left( x = || m \|| \right)$ / $\pi \left[ x = \phi \left( m \right) \right]$",
                                  fontsize=15, labelpad=10)
            plt.tick_params(axis="both", which="major", labelsize=10, pad=10)

            axis[0, 0].hist(cartesian_magnetisation_sample[0, no_of_equilibration_sweeps:],
                            bins=number_of_histogram_bins, density=True)
            axis[0, 1].hist(cartesian_magnetisation_sample[1, no_of_equilibration_sweeps:],
                            bins=number_of_histogram_bins, density=True)
            axis[1, 0].hist(np.abs(cartesian_magnetisation_sample[0, no_of_equilibration_sweeps:]),
                            bins=number_of_histogram_bins, density=True)
            axis[1, 1].hist(np.abs(cartesian_magnetisation_sample[1, no_of_equilibration_sweeps:]),
                            bins=number_of_histogram_bins, density=True)
            axis[2, 0].hist(magnetisation_norm_sample[no_of_equilibration_sweeps:], bins=number_of_histogram_bins,
                            density=True)
            axis[2, 1].hist(magnetisation_phase_sample[no_of_equilibration_sweeps:], bins=number_of_histogram_bins,
                            density=True)
            plt.savefig(f"{output_directory}/magnetisation_histograms_temp_eq_{temperature:.2f}_"
                        f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}_"
                        f"job_{job_index + 1}.pdf", bbox_inches="tight")
            plt.clf()

        temperature -= magnitude_of_temperature_increments
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        raise Exception("InterfaceError: One positional argument required - give the configuration-file location.  In "
                        "addition, you may provide number_of_histogram_bins (default value is 1000) in the second "
                        "position.")
    if len(sys.argv) == 2:
        print("One positional argument provided.  In addition, you may provide number_of_histogram_bins (default value "
              "is 1000) in the second position.")
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        print("Two positional arguments provided.  The second must be number_of_histogram_bins.")
        main(sys.argv[1], int(sys.argv[2]))
