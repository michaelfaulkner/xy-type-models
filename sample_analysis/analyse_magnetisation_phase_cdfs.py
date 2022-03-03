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
    (algorithm_name, output_directory, no_of_sites, no_of_equilibration_sweeps, no_of_observations, temperatures,
     use_external_global_moves, external_global_moves_string, no_of_jobs,
     max_no_of_cpus) = run_script.get_config_data(config_file)
    if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
        print("ConfigurationError: The configuration file corresponds to a Maggs-electrolyte model but this script "
              "requires the XY of HXY model.")
        raise SystemExit
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
    colors = ["red", "black"]
    linestyles = ["solid", "dotted", "dashed", "dashdot", (0, (1, 1)), (0, (5, 10)), (0, (5, 1)), (0, (3, 1, 1, 1))]
    inset_axis = plt.axes([0.135, 0.665, 0.2, 0.2])
    inset_axis.yaxis.set_label_position("right")
    inset_axis.yaxis.tick_right()
    inset_axis.tick_params(which='major', width=2, labelsize=7.5)
    inset_axis.set_xlabel(r"$1 / (\beta J)$", fontsize=7.5)
    inset_axis.set_ylabel(r"$\omega_n^2$", fontsize=7.5)
    [inset_axis.spines[spine].set_linewidth(2) for spine in ["top", "bottom", "left", "right"]]

    try:
        with open(f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}_"
                  f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites_"
                  f"{no_of_observations}_obs.tsv", "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            cvms, cvm_errors = output_file_sans_header[1], output_file_sans_header[2]
        for temperature_index, temperature in setup_scripts.reverse_enumerate(temperatures):
            if temperature_index == 0 or temperature_index == len(temperatures) - 1:
                if no_of_jobs == 1:
                    axis.plot(*np.load(f"{output_directory}/magnetisation_phase_cdf_{algorithm_name.replace('-', '_')}_"
                                       f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x"
                                       f"{int(no_of_sites ** 0.5)}_sites_{no_of_observations}_obs_temp_eq_"
                                       f"{temperature:.4f}.npy"), color=colors[min(temperature_index, 1)],
                              label=f"temperature = {temperature:.4f}")
                else:
                    for job_index in range(min(no_of_jobs, no_of_plotted_cdfs)):
                        output_file = (f"{output_directory}/magnetisation_phase_cdf_{algorithm_name.replace('-', '_')}_"
                                       f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x"
                                       f"{int(no_of_sites ** 0.5)}_sites_{no_of_observations}_obs_temp_eq_"
                                       f"{temperature:.4f}_job_{job_index}.npy")
                        if job_index == 0:
                            axis.plot(*np.load(output_file), color=colors[min(temperature_index, 1)],
                                      linestyle=linestyles[0], label=fr"$1 / (\beta J)$ = {temperature:.4f}")
                            """use alternative CDF calculation in 
                                markov_chain_diagnostics.get_cumulative_distribution() to use fill_between() in the 
                                following line"""
                            """if temperature_index != 0:
                                axis.fill_between(cdf[0], cdf[1], np.arange(0.0, 1.0, 1.0 / float(len(cdf[0]))), 
                                                  color='green')"""
                        else:
                            axis.plot(*np.load(output_file), color=colors[min(temperature_index, 1)],
                                      linestyle=linestyles[job_index])
    except IOError:
        cvm_file = open(f"{output_directory}/magnetisation_phase_cramervonmises_{algorithm_name.replace('-', '_')}_"
                        f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites_"
                        f"{no_of_observations}_obs.tsv", "w")
        cvm_file.write("# temperature".ljust(30) + "omega_n^2".ljust(30) + "omega_n^2 error" + "\n")
        cvms, cvm_errors = [], []
        start_time = time.time()
        for temperature_index, temperature in setup_scripts.reverse_enumerate(temperatures):
            print(f"Temperature = {temperature:.4f}")
            if no_of_jobs == 1:
                cdf_of_magnetisation_phase, cvm = get_cdf_and_cramervonmises_of_magnetisation_phase(
                    output_directory, temperature, temperature_index, no_of_sites, no_of_equilibration_sweeps)
                cvm_error = 1.0  # dummy error
                if temperature_index == 0 or temperature_index == len(temperatures) - 1:
                    axis.plot(*cdf_of_magnetisation_phase, color=colors[min(temperature_index, 1)],
                              label=f"temperature = {temperature:.4f}")
                    np.save(f"{output_directory}/magnetisation_phase_cdf_{algorithm_name.replace('-', '_')}_"
                            f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites_"
                            f"{no_of_observations}_obs_temp_eq_{temperature:.4f}.npy", cdf_of_magnetisation_phase)
            else:
                cdf_and_cvm = np.array(pool.starmap(get_cdf_and_cramervonmises_of_magnetisation_phase, [
                    (f"{output_directory}/job_{job_index}", temperature, temperature_index, no_of_sites,
                     no_of_equilibration_sweeps) for job_index in range(no_of_jobs)]), dtype=object).transpose()
                cdfs_of_magnetisation_phase, cvm, cvm_error = cdf_and_cvm[0], np.mean(cdf_and_cvm[1]), (
                        np.std(cdf_and_cvm[1]) / len(cdf_and_cvm[1]) ** 0.5)
                if temperature_index == 0 or temperature_index == len(temperatures) - 1:
                    for job_index, cdf in enumerate(cdfs_of_magnetisation_phase[:min(no_of_jobs, no_of_plotted_cdfs)]):
                        if job_index == 0:
                            axis.plot(*cdf, color=colors[min(temperature_index, 1)], linestyle=linestyles[0],
                                      label=fr"$1 / (\beta J)$ = {temperature:.4f}")
                        else:
                            axis.plot(*cdf, color=colors[min(temperature_index, 1)], linestyle=linestyles[job_index])
                        np.save(f"{output_directory}/magnetisation_phase_cdf_{algorithm_name.replace('-', '_')}_"
                                f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_"
                                f"sites_{no_of_observations}_obs_temp_eq_{temperature:.4f}_job_{job_index}.npy", cdf)
            cvms.append(cvm)
            cvm_errors.append(cvm_error)
            cvm_file.write(f"{temperature:.14e}".ljust(30) + f"{cvm:.14e}".ljust(30) + f"{cvm_error:.14e}" + "\n")

        cvm_file.close()
        if no_of_jobs > 1:
            pool.close()
        print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    if no_of_jobs == 1:
        inset_axis.errorbar(list(reversed(temperatures)), cvms, marker=".", markersize=10, color="k", linestyle="None")
    else:
        inset_axis.errorbar(list(reversed(temperatures)), cvms, cvm_errors, marker=".", markersize=10, color="k",
                            linestyle="None")

    handles, labels = axis.get_legend_handles_labels()
    legend = axis.legend(reversed(handles), reversed(labels), title='colour code', loc="lower right", fontsize=10)
    legend.get_frame().set_edgecolor("k")
    legend.get_frame().set_lw(1.5)
    figure.savefig(f"{output_directory}/magnetisation_phase_CDFs_w_cramervonmises_{algorithm_name.replace('-', '_')}_"
                   f"{external_global_moves_string}_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites_"
                   f"{no_of_observations}_obs.pdf", bbox_inches="tight")


def get_cdf_and_cramervonmises_of_magnetisation_phase(sample_directory, temperature, temperature_index, no_of_sites,
                                                      no_of_equilibration_sweeps):
    r"""
    Returns (for the sample of the magnetisation phase) the empirical CDF F_n(x) and normalised Cramér-von Mises
    integral \omega_n^2 := \int (F_n(x) - F(x))^2 dF(x) / n, where n is the sample size and F(x) is the CDF of the
    continuous uniform distribution Uniform(-pi, pi).

    The magnetisation phase phi(x; temperature, no_of_sites), where
    m(x; temperature, no_of_sites) = (|| m(x; temperature, no_of_sites) ||, phi(x; temperature, no_of_sites))^t in
    radial coordinates and m(x; temperature, no_of_sites) = sum_i [cos(x_i), sin(x_i)]^t / no_of_sites is the Cartesian
    magnetisation, with x_i the position of particle i at the time of observation.

    NB, RuntimeWarnings are often thrown by lines 175, 178, 179 and 188 of scipy/stats/_hypotests.py when calculating
    the p-value of scipy.stats.cramervonmises() for highly nonergodic magnetisation-phase samples.  They do not affect
    the calculation of the Cramér-von Mises integral itself.

    Parameters
    ----------
    sample_directory : str
        The location of the directory containing the sample and Metropolis acceptance rate.
    temperature : float
        The sampling temperature.
    temperature_index : int
        The index of the current sampling temperature within the configuration file.
    no_of_sites : int
        The number of lattice sites.
    no_of_equilibration_sweeps : int
        The number of discarded equilibration observations.

    Returns
    -------
    numpy.ndarray
        the empirical CDF and normalised Cramér-von Mises integral.  A one-dimensional Python list of length 2.  The
        first element is a two-dimensional numpy array (of floats) of size (2, n) representing the empirical CDF.  The
        second element is a float estimating the normalised Cramér-von Mises integral.
    """
    magnetisation_phase = sample_getter.get_magnetisation_phase(sample_directory, temperature, temperature_index,
                                                                no_of_sites)[no_of_equilibration_sweeps:]
    return [np.array(markov_chain_diagnostics.get_cumulative_distribution(magnetisation_phase)),
            stats.cramervonmises(magnetisation_phase, cdf="uniform", args=(-math.pi, 2.0 * math.pi)).statistic /
            len(magnetisation_phase)]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise Exception("InterfaceError: One positional argument required - give the configuration-file location.")
    else:
        print("One positional argument provided.  It must be the location of the configuration file.")
        main(sys.argv[1])
