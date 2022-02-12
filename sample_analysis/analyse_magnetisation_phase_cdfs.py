from scipy import stats
import importlib
import math
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


def main(config_file, no_of_plotted_cdfs=8):
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
    axis.set_xlim([-math.pi, math.pi])
    axis.set_xticks(np.arange(-math.pi, math.pi + 0.5 * math.pi / 2, step=(0.5 * math.pi)))
    axis.set_xticklabels([r"$-\pi$", r"$-\pi / 2$", r"$0$", r"$\pi / 2$", r"$\pi$"])
    axis.set_xlabel(r"$x$", fontsize=20, labelpad=8)
    axis.set_ylim(0.0, 1.0)
    axis.tick_params(which='major', width=2, length=7, labelsize=18, pad=10)
    axis.set_ylabel(r"$\mathbb{P}_n \left( \phi_m < x \right)$", fontsize=20, labelpad=-30)
    axis.yaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
    axis.yaxis.set_major_formatter('{x:.1f}')
    [axis.spines[spine].set_linewidth(2) for spine in ["top", "bottom", "left", "right"]]
    colors = ["black", "red"]
    linestyles = ["solid", "dotted", "dashed", "dashdot", (0, (1, 1)), (0, (5, 10)), (0, (5, 1)), (0, (3, 1, 1, 1))]
    inset_axis = plt.axes([0.15, 0.6125, 0.25, 0.25])
    inset_axis.yaxis.set_label_position("right")
    inset_axis.yaxis.tick_right()
    inset_axis.set_xlabel(r"$1 / (\beta J)$")
    inset_axis.set_ylabel(r"$\omega_{\rm{CvM}}$")
    [inset_axis.spines[spine].set_linewidth(2) for spine in ["top", "bottom", "left", "right"]]

    try:
        with open(f"{output_directory}/cramervonmises_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_"
                  f"{algorithm_name.replace('-', '_')}.tsv", "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            temperatures, cvms, cvm_errors = (output_file_sans_header[0], output_file_sans_header[1],
                                              output_file_sans_header[2])
        for temperature_index in range(min(no_of_temperature_increments + 1, 2)):
            if no_of_jobs == 1:
                with open(f"{output_directory}/magnetisation_phase_cdf_temp_eq_{temperature:.2f}_"
                          f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_"
                          f"{algorithm_name.replace('-', '_')}.npy", "rb") as output_file:
                    axis.plot(*np.load(output_file), color=colors[min(temperature_index, 1)],
                              label=f"temperature = {temperature:.2f}")
            else:
                for job_index in range(min(no_of_jobs, no_of_plotted_cdfs)):
                    with open(f"{output_directory}/magnetisation_phase_cdf_temp_eq_{temperature:.2f}_"
                              f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}_"
                              f"job_{job_index + 1}.npy", "rb") as output_file:
                        if job_index == 0:
                            axis.plot(*np.load(output_file), color=colors[min(temperature_index, 1)],
                                      linestyle=linestyles[0], label=fr"$1 / (\beta J)$ = {temperature:.2f}")
                            """use alternative CDF calculation in markov_chain_diagnostics.get_cumulative_distribution()
                                to use fill_between() in the following line"""
                            """if temperature_index != 0:
                                axis.fill_between(cdf[0], cdf[1], np.arange(0.0, 1.0, 1.0 / float(len(cdf[0]))), 
                                                  color='green')"""
                        else:
                            axis.plot(*np.load(output_file), color=colors[min(temperature_index, 1)],
                                      linestyle=linestyles[job_index])
            temperature = initial_temperature
    except IOError:
        cvm_file = open(f"{output_directory}/cramervonmises_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_"
                        f"{algorithm_name.replace('-', '_')}.tsv", "w")
        cvm_file.write("# temperature".ljust(30) + "Cramer-von Mises value".ljust(30) + "Cramer-von Mises error" + "\n")
        temperatures, cvms, cvm_errors = [], [], []
        start_time = time.time()
        for temperature_index in range(no_of_temperature_increments + 1):
            print(f"Temperature = {temperature:.2f}")
            if no_of_jobs == 1:
                cdf_of_magnetisation_phase, cvm_squared = get_cdf_and_cramervonmises_sq_of_magnetisation_phase(
                    output_directory, temperature, no_of_sites, no_of_equilibration_sweeps)
                cvm, cvm_error = cvm_squared ** 0.5, 1.0  # dummy error
                if temperature_index == 0 or temperature_index == no_of_temperature_increments:
                    axis.plot(*cdf_of_magnetisation_phase, color=colors[min(temperature_index, 1)],
                              label=f"temperature = {temperature:.2f}")
                    with open(f"{output_directory}/magnetisation_phase_cdf_temp_eq_{temperature:.2f}_"
                              f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_"
                              f"{algorithm_name.replace('-', '_')}.npy", "wb") as output_file:
                        np.save(output_file, cdf_of_magnetisation_phase)
            else:
                cdf_and_cvm_sq = np.array(pool.starmap(get_cdf_and_cramervonmises_sq_of_magnetisation_phase, [
                    (f"{output_directory}/job_{job_index + 1}", temperature, no_of_sites, no_of_equilibration_sweeps)
                    for job_index in range(no_of_jobs)]), dtype=object).transpose()
                cdfs_of_magnetisation_phase, cvm, cvm_error = cdf_and_cvm_sq[0], np.mean(cdf_and_cvm_sq[1]) ** 0.5, (
                        np.std(cdf_and_cvm_sq[1]) / len(cdf_and_cvm_sq[1]) ** 0.5) ** 0.5
                if temperature_index == 0 or temperature_index == no_of_temperature_increments:
                    for job_index, cdf in enumerate(cdfs_of_magnetisation_phase[:min(no_of_jobs, no_of_plotted_cdfs)]):
                        if job_index == 0:
                            axis.plot(*cdf, color=colors[min(temperature_index, 1)], linestyle=linestyles[0],
                                      label=fr"$1 / (\beta J)$ = {temperature:.2f}")
                        else:
                            axis.plot(*cdf, color=colors[min(temperature_index, 1)], linestyle=linestyles[job_index])
                        with open(f"{output_directory}/magnetisation_phase_cdf_temp_eq_{temperature:.2f}_"
                                  f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_"
                                  f"{algorithm_name.replace('-', '_')}_job_{job_index + 1}.npy", "wb") as output_file:
                            np.save(output_file, cdf)
            temperatures.append(temperature)
            cvms.append(cvm)
            cvm_errors.append(cvm_error)
            cvm_file.write(f"{temperature:.14e}".ljust(30) + f"{cvm:.14e}".ljust(30) + f"{cvm_error:.14e}" + "\n")
            temperature -= magnitude_of_temperature_increments

        cvm_file.close()
        if no_of_jobs > 1:
            pool.close()
        print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    if no_of_jobs == 1:
        inset_axis.plot(temperatures, cvms, marker=".", markersize=10, color="k", linestyle="None")
    else:
        inset_axis.errorbar(temperatures, cvms, cvm_errors, marker=".", markersize=10, color="k", linestyle="None")

    handles, labels = axis.get_legend_handles_labels()
    legend = axis.legend(reversed(handles), reversed(labels), title='Colour code', loc="lower right", fontsize=10)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(1.5)

    if use_external_global_moves:
        figure.savefig(f"{output_directory}/magnetisation_phase_CDFs_w_twists_{int(no_of_sites ** 0.5)}x"
                       f"{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf", bbox_inches="tight")
    else:
        figure.savefig(f"{output_directory}/magnetisation_phase_CDFs_sans_twists_{int(no_of_sites ** 0.5)}x"
                       f"{int(no_of_sites ** 0.5)}_{algorithm_name.replace('-', '_')}.pdf", bbox_inches="tight")


def get_cdf_and_cramervonmises_sq_of_magnetisation_phase(sample_directory, temperature, no_of_sites,
                                                         no_of_equilibration_sweeps):
    magnetisation_phase = sample_getter.get_magnetisation_phase(sample_directory, temperature,
                                                                no_of_sites)[no_of_equilibration_sweeps:]
    return [np.array(markov_chain_diagnostics.get_cumulative_distribution(magnetisation_phase)),
            stats.cramervonmises(magnetisation_phase, cdf="uniform", args=(-math.pi, 2.0 * math.pi)).statistic]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise Exception("InterfaceError: One positional argument required - give the configuration-file location.")
    else:
        print("One positional argument provided.  It must be the location of the configuration file.")
        main(sys.argv[1])
