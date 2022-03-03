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


def main(config_file, observable_string, max_physical_time=100.0):
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, _, temperatures, _,
     external_global_moves_string, no_of_jobs, max_no_of_cpus) = run_script.get_config_data(config_file)
    setup_scripts.check_for_observable_error(algorithm_name, observable_string)
    get_sample_method = getattr(sample_getter, "get_" + observable_string)
    linestyles = ["solid", "dotted", "dashed", "dashdot", (0, (1, 1)), (0, (5, 10)), (0, (5, 1)), (0, (3, 1, 1, 1))]
    sample_is_one_dimensional = setup_scripts.get_sample_is_one_dimensional(observable_string)

    start_time = time.time()
    for temperature_index, temperature in enumerate(temperatures):
        print(f"Temperature = {temperature:.4f}")
        if no_of_jobs == 1:
            physical_time_step = sample_getter.get_physical_time_step(algorithm_name, f"{output_directory}",
                                                                      temperature_index)
        else:
            physical_time_step = sum([sample_getter.get_physical_time_step(
                algorithm_name, f"{output_directory}/job_{job_index}", temperature_index)
                for job_index in range(no_of_jobs)]) / no_of_jobs
        plt.figure(0)
        plt.xlabel(r"total simulation time, $t_{\rm{tot}}$ ($s$)", fontsize=15, labelpad=10)
        plt.ylabel(r"$\bar{x}(t_{\rm{tot}})$ (red) / $\sigma_x^2(t_{\rm{tot}})$ (blue)", fontsize=10, labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=14, pad=10)

        for job_index in range(no_of_jobs):
            try:
                # first, attempt to open file containing the sample, running mean and running variance...
                sample, running_mean, running_variance = np.load(
                    f"{output_directory}/{observable_string}_sample_and_running_mean_and_variance_"
                    f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{int(no_of_sites ** 0.5)}x"
                    f"{int(no_of_sites ** 0.5)}_sites_temp_eq_{temperature:.4f}_job_{job_index}.npy")
            except IOError:
                # ...then compute the sample, running mean and running variance if the file does not exist
                if no_of_jobs > 1:
                    sample_directory = f"{output_directory}/job_{job_index}"
                    sample = get_sample_method(sample_directory, temperature, temperature_index, no_of_sites)[
                             no_of_equilibration_sweeps:]
                else:
                    sample = get_sample_method(output_directory, temperature, temperature_index, no_of_sites)[
                             no_of_equilibration_sweeps:]
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
            reduced_sample = sample[:int(max_total_physical_time / physical_time_step) + 1]
            reduced_running_mean = running_mean[:int(max_total_physical_time / physical_time_step) + 1]
            reduced_running_variance = running_variance[:int(max_total_physical_time / physical_time_step) + 1]
            '''if temperature_index == 0:
                reduced_sample = sample
                reduced_running_mean = running_mean
                reduced_running_variance = running_variance
            else:
                reduced_sample_length = int(len(sample))
                reduced_sample = sample[:reduced_sample_length + 1]
                reduced_running_mean = running_mean[:reduced_sample_length + 1]
                reduced_running_variance = running_variance[:reduced_sample_length + 1]'''

            plt.figure(0)
            plt.plot(np.arange(len(reduced_running_mean)) * physical_time_step, reduced_running_mean, color="r",
                     linewidth=1, linestyle=linestyles[job_index])
            plt.plot(np.arange(len(reduced_running_variance)) * physical_time_step, reduced_running_variance, color="b",
                     linewidth=1, linestyle=linestyles[job_index])

            plt.figure(job_index + 1)
            plt.xlabel(r"time, $t$ ($s$)", fontsize=15, labelpad=10)
            plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
            if sample_is_one_dimensional:
                plt.ylabel(r"$ x \left( t \right) - \bar{x}$", fontsize=15, labelpad=10)
            else:
                plt.ylabel(r"$ x_1 \left( t \right) - \bar{x_1}$", fontsize=15, labelpad=10)
            plt.plot(np.arange(len(reduced_sample)) * physical_time_step, reduced_sample, color="k", linewidth=1,
                     linestyle="-")
            plt.plot(np.arange(len(reduced_running_mean)) * physical_time_step, reduced_running_mean, color="r",
                     linewidth=1, linestyle="-")
            plt.plot(np.arange(len(reduced_running_variance)) * physical_time_step, reduced_running_variance, color="b",
                     linewidth=1, linestyle="-")
            plt.tight_layout()
            plt.savefig(f"{output_directory}/{observable_string}_vs_time_{algorithm_name.replace('-', '_')}_"
                        f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites_"
                        f"temp_eq_{temperature:.4f}_job_{job_index}.pdf", bbox_inches="tight")
            # plt.clf()

        plt.figure(0)
        plt.tight_layout()
        plt.savefig(f"{output_directory}/{observable_string}_running_mean_and_variance_vs_time_"
                    f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{int(no_of_sites ** 0.5)}x"
                    f"{int(no_of_sites ** 0.5)}_sites_temp_eq_{temperature:.4f}.pdf", bbox_inches="tight")
        plt.clf()
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose trace plot and histogram you wish to estimate (in the first"
                        " and second positions, respectively).  In addition, you may provide max_physical_time (default"
                        " value is 100.0) in the third and fourth position.")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  The first / second must be the location of the configuration file / "
              "the string of the observable whose trace plot and histogram you wish to estimate.  In addition, you may "
              "provide max_physical_time (default value is 100.0) in the third position.")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The first / second / third must be the location of the "
              "configuration file / the string of the observable whose trace plot and histogram you wish to estimate / "
              "max_physical_time.")
        main(sys.argv[1], sys.argv[2], float(sys.argv[3]))

