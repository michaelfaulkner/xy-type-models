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


def main(observable_string="rotated_magnetisation_phase"):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    config_file = "config_files/running_mean_figure.txt"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, _, temperatures, _,
     external_global_moves_string, no_of_jobs, max_no_of_cpus) = run_script.get_config_data(config_file)
    setup_scripts.check_for_observable_error(algorithm_name, observable_string)
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
    linestyles = ["solid", "dotted", "dashed", "dashdot", (0, (1, 1)), (0, (5, 10)), (0, (5, 1)), (0, (3, 1, 1, 1))]
    sample_is_one_dimensional = setup_scripts.get_sample_is_one_dimensional(observable_string)

    figure, axes = plt.subplots(2, 1, figsize=(7.5, 7.1))
    axes[1].set_xlabel(r"$\tau$", fontsize=20, labelpad=0)
    axes[0].set_ylabel(r"$\bar{\phi}_m(\tau)$ / $s_{\phi_m}^2$", fontsize=20, labelpad=-10)
    axes[1].set_ylabel(r"$\bar{\phi}_m(\tau)$", fontsize=20, labelpad=-10)
    [axis.tick_params(which='major', width=3, length=7, labelsize=18) for axis in axes]
    [axis.locator_params(axis='x', nbins=4) for axis in axes]
    [axis.locator_params(axis='y', nbins=8) for axis in axes]
    axes[0].xaxis.set_tick_params(labelbottom=False)
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes]
    colors = ["black", "red"]

    start_time = time.time()
    for temperature_index, temperature in enumerate(temperatures):
        print(f"Temperature = {temperature:.4f}")
        try:
            physical_time_step = np.loadtxt(f"{output_directory}/physical_time_step_{algorithm_name.replace('-', '_')}_"
                                            f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x"
                                            f"{int(no_of_sites ** 0.5)}_sites_temp_eq_{temperature:.4f}.csv",
                                            dtype=float, delimiter=",")
        except IOError:
            physical_time_step = sum([sample_getter.get_physical_time_step(
                    algorithm_name, f"{output_directory}/job_{job_index}", temperature_index) for job_index in
                range(no_of_jobs)]) / no_of_jobs
            np.savetxt(f"{output_directory}/physical_time_step_{algorithm_name.replace('-', '_')}_"
                       f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites_"
                       f"temp_eq_{temperature:.4f}.csv", np.atleast_1d(physical_time_step), delimiter=",")

        for job_index in range(no_of_jobs):
            try:
                # first, attempt to open file containing the sample, running mean and running variance...
                sample, running_mean, running_variance = np.load(
                    f"{output_directory}/{observable_string}_sample_and_running_mean_and_variance_"
                    f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{int(no_of_sites ** 0.5)}x"
                    f"{int(no_of_sites ** 0.5)}_sites_temp_eq_{temperature:.4f}_job_{job_index}.npy")
            except IOError:
                # ...then compute the sample, running mean and running variance if the file does not exist
                sample = get_sample_method(f"{output_directory}/job_{job_index}", temperature, temperature_index,
                                           no_of_sites)[no_of_equilibration_sweeps:]
                if not sample_is_one_dimensional:
                    sample = sample.transpose()[0]
                running_mean = np.array([np.mean(sample[:index + 1]) for index in range(len(sample))])
                running_variance = np.array([np.var(sample[:index + 1]) for index in range(len(sample))])
                np.save(f"{output_directory}/{observable_string}_sample_and_running_mean_and_variance_"
                        f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{int(no_of_sites ** 0.5)}x"
                        f"{int(no_of_sites ** 0.5)}_sites_temp_eq_{temperature:.4f}_job_{job_index}.npy",
                        np.array([sample, running_mean, running_variance]))

            if temperature_index == 0 and job_index == 0:
                max_total_physical_time = physical_time_step * (len(sample) - 1)
            reduced_running_variance = running_variance[:int(max_total_physical_time / physical_time_step) + 1]
            reduced_running_mean = (running_mean[:int(max_total_physical_time / physical_time_step) + 1] /
                                    reduced_running_variance[-1])

            if job_index == 0:
                axes[0].plot(np.arange(len(reduced_running_mean)) * physical_time_step, reduced_running_mean,
                             color=colors[temperature_index], linewidth=2, linestyle=linestyles[job_index],
                             label=fr"$1 / (\beta J)$ = {temperature:.2f}")
                axes[1].plot(np.arange(len(reduced_running_mean)) * physical_time_step, reduced_running_mean,
                             color=colors[temperature_index], linewidth=2, linestyle=linestyles[job_index],
                             label=fr"$1 / (\beta J)$ = {temperature:.2f}")
            else:
                axes[0].plot(np.arange(len(reduced_running_mean)) * physical_time_step, reduced_running_mean,
                             color=colors[temperature_index], linewidth=2, linestyle=linestyles[job_index])
                axes[1].plot(np.arange(len(reduced_running_mean)) * physical_time_step, reduced_running_mean,
                             color=colors[temperature_index], linewidth=2, linestyle=linestyles[job_index])

    figure.tight_layout()
    legend = axes[0].legend(loc="upper left", fontsize=12)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(3)
    figure.savefig(f"{output_directory}/{observable_string}_running_mean_vs_time_{algorithm_name.replace('-', '_')}_"
                   f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites.pdf",
                   bbox_inches="tight")
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("InterfaceError: No positional arguments required.  You may provide the string of the "
                        "observable whose running sample mean you wish to plot (default is "
                        "'rotated_magnetisation_phase').")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be the the string of the observable whose running sample "
              "mean you wish to plot (default is 'rotated_magnetisation_phase'.")
        main(sys.argv[1])
    else:
        main()
