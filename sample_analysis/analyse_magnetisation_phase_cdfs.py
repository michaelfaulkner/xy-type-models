import importlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import sys
import time

# import additional modules; have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")
setup_scripts = importlib.import_module("setup_scripts")
sample_getter = importlib.import_module("sample_getter")
run_script = importlib.import_module("run")


def main(config_file, no_of_histogram_bins=100):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, initial_temperature, final_temperature,
     no_of_temperature_increments, use_external_global_moves, no_of_jobs, max_no_of_cpus) = run_script.get_config_data(
        config_file)
    if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
        print("ConfigurationError: The configuration file corresponds to a Maggs-electrolyte model but this script "
              "requires the XY of HXY model.")
        raise SystemExit
    (temperature, magnitude_of_temperature_increments) = setup_scripts.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)
    pool = setup_scripts.setup_pool(no_of_jobs, max_no_of_cpus)
    figure, axis = plt.subplots(1)
    axis.tick_params(which='both', width=2)
    axis.tick_params(which='major', length=7, labelsize=18, pad=10)
    axis.tick_params(which='minor', length=4)
    axis.set_ylim(0.0, 1.0)
    axis.tick_params(which='major', width=2, length=7, labelsize=18, pad=10)
    axis.set_ylabel(r"$\mathbb{P} \left( \phi_m < x \right)$", fontsize=20, labelpad=-30)
    axis.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    axis.yaxis.set_major_formatter('{x:.1f}')
    colors = plt.cm.rainbow(np.linspace(0, 1, no_of_temperature_increments + 1))
    linestyles = ["-", "--", "-.", ":"]

    start_time = time.time()
    for temperature_index in range(no_of_temperature_increments):
        print(f"Temperature = {temperature:.2f}")
        if no_of_jobs == 1:
            cdfs_of_magnetisation_phase = np.atleast_2d(
                get_cdf_of_magnetisation_phase(output_directory, temperature, no_of_sites, no_of_equilibration_sweeps))
        else:
            cdfs_of_magnetisation_phase = np.array(pool.starmap(get_cdf_of_magnetisation_phase, [
                (f"{output_directory}/job_{job_index + 1}", temperature, no_of_sites, no_of_equilibration_sweeps)
                for job_index in range(no_of_jobs)]))
        if temperature_index == 0 or temperature_index == no_of_temperature_increments - 1:
            for index, cdf in enumerate(cdfs_of_magnetisation_phase[:min(no_of_jobs, 4)]):
                axis.plot(*cdf, color=colors[temperature_index], linestyle=linestyles[index])
        temperature += magnitude_of_temperature_increments
    pool.close()

    if use_external_global_moves:
        figure.savefig(f"{output_directory}/magnetisation_phase_CDFs_w_twists_temp_eq_{temperature:.2f}_"
                       f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}"
                       f".pdf", bbox_inches="tight")
    else:
        figure.savefig(f"{output_directory}/magnetisation_phase_CDFs_w_twists_temp_eq_{temperature:.2f}_"
                       f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}"
                       f".pdf", bbox_inches="tight")
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


def get_cdf_of_magnetisation_phase(sample_directory, temperature, no_of_sites, no_of_equilibration_sweeps):
    return markov_chain_diagnostics.get_cumulative_distribution(
        sample_getter.get_magnetisation_phase(sample_directory, temperature, no_of_sites)[no_of_equilibration_sweeps:])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise Exception("InterfaceError: One positional argument required - give the configuration-file location.")
    else:
        print("One positional argument provided.  It must be the location of the configuration file.")
        main(sys.argv[1])
