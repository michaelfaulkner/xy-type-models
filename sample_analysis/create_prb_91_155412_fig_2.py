import importlib
import math
import matplotlib.pyplot as plt
import multiprocessing as mp
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
markov_chain_diagnostics = importlib.import_module("markov_chain_diagnostics")
run_script = importlib.import_module("run")


def main(electrolyte_version):
    if electrolyte_version:
        config_file_local_moves = "config_files/prb_91_155412_fig_2/local_moves_only.txt"
        config_file_all_moves = "config_files/prb_91_155412_fig_2/all_moves.txt"
    else:
        config_file_local_moves = "config_files/prb_91_155412_fig_2_hxy/local_moves_only.txt"
        config_file_all_moves = "config_files/prb_91_155412_fig_2_hxy/all_moves.txt"
    (algorithm_name, output_directory_local_moves, no_of_sites, no_of_equilibration_sweeps, initial_temperature,
     final_temperature, no_of_temperature_increments, _, no_of_jobs, max_no_of_cpus) = run_script.get_config_data(
        config_file_local_moves)
    output_directory = output_directory_local_moves.replace("/local_moves_only", "")
    output_directory_all_moves = run_script.get_config_data(config_file_all_moves)[1]
    (temperature, magnitude_of_temperature_increments) = setup_scripts.get_temperature_and_magnitude_of_increments(
        initial_temperature, final_temperature, no_of_temperature_increments)

    try:
        with open(f"{output_directory}/prb_91_155412_fig_2_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_"
                  f"{algorithm_name.replace('-', '_')}.tsv", "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            temperatures = output_file_sans_header[0]
            chi_w_ratios = output_file_sans_header[1]
    except IOError:
        output_file = open(f"{output_directory}/prb_91_155412_fig_2_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_"
                           f"{algorithm_name.replace('-', '_')}.tsv", "w")

        if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
            output_file.write("# temperature".ljust(15) + "chi_w ratio".ljust(30) + "chi_w ratio error".ljust(30) +
                              "final width of proposal interval (local only)".ljust(50) +
                              "rotational acceptance rate (local only)".ljust(50) +
                              "charge-hop acceptance rate (local only)".ljust(50) +
                              "final width of proposal interval (all moves)".ljust(50) +
                              "rotational acceptance rate (all moves)".ljust(50) +
                              "charge-hop acceptance rate (all moves)".ljust(50) +
                              "acceptance rate (external global moves)" + "\n")
        else:
            output_file.write("# temperature".ljust(15) + "chi_w ratio".ljust(30) + "chi_w ratio error".ljust(30) +
                              "final width of proposal interval (local only)".ljust(50) +
                              "rotational acceptance rate (local only)".ljust(50) +
                              "final width of proposal interval (all moves)".ljust(50) +
                              "rotational acceptance rate (all moves)".ljust(50) +
                              "acceptance rate (external global moves)" + "\n")

        no_of_cpus = mp.cpu_count()
        pool = mp.Pool(no_of_cpus)
        temperatures = []
        chi_w_ratios = []
        start_time = time.time()

        for i in range(no_of_temperature_increments + 1):
            print(f"Temperature = {temperature:.2f}")

            acceptance_rates_local_moves = np.mean(
                [sample_getter.get_acceptance_rates(output_directory_local_moves + "/job_" + str(job_number + 1),
                                                    temperature) for job_number in range(no_of_jobs)], axis=0)
            acceptance_rates_all_moves = np.mean(
                [sample_getter.get_acceptance_rates(output_directory_all_moves + "/job_" + str(job_number + 1),
                                                    temperature) for job_number in range(no_of_jobs)], axis=0)

            sample_means_and_errors_local_moves = np.transpose(
                np.array(pool.starmap(markov_chain_diagnostics.get_sample_mean_and_error, [[
                    sample_getter.get_topological_susceptibility(
                        output_directory_local_moves + "/job_" + str(job_number + 1), temperature, no_of_sites)[
                        no_of_equilibration_sweeps:]] for job_number in range(no_of_jobs)])))
            sample_means_and_errors_all_moves = np.transpose(
                np.array(pool.starmap(markov_chain_diagnostics.get_sample_mean_and_error, [[
                    sample_getter.get_topological_susceptibility(
                        output_directory_all_moves + "/job_" + str(job_number + 1), temperature, no_of_sites)[
                        no_of_equilibration_sweeps:]] for job_number in range(no_of_jobs)])))

            sample_mean_local_moves = np.mean(sample_means_and_errors_local_moves[0])
            sample_error_local_moves = np.linalg.norm(sample_means_and_errors_local_moves[1])
            sample_mean_all_moves = np.mean(sample_means_and_errors_all_moves[0])
            sample_error_all_moves = np.linalg.norm(sample_means_and_errors_all_moves[1])
            chi_w_ratio = sample_mean_local_moves / sample_mean_all_moves
            chi_w_ratio_error = math.sqrt((sample_error_local_moves / sample_mean_all_moves) ** 2 +
                                          (sample_mean_local_moves * sample_error_all_moves /
                                           sample_mean_all_moves ** 2) ** 2)

            if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
                output_file.write(f"{temperature:.2f}".ljust(15) + f"{chi_w_ratio:.14e}".ljust(30) +
                                  f"{chi_w_ratio_error:.14e}".ljust(30) +
                                  f"{acceptance_rates_local_moves[0]:.14e}".ljust(50) +
                                  f"{acceptance_rates_local_moves[1]:.14e}".ljust(50) +
                                  f"{acceptance_rates_local_moves[2]:.14e}".ljust(50) +
                                  f"{acceptance_rates_all_moves[0]:.14e}".ljust(50) +
                                  f"{acceptance_rates_all_moves[1]:.14e}".ljust(50) +
                                  f"{acceptance_rates_all_moves[2]:.14e}".ljust(50) +
                                  f"{acceptance_rates_all_moves[3]:.14e}" + "\n")
            else:
                output_file.write(f"{temperature:.2f}".ljust(15) + f"{chi_w_ratio:.14e}".ljust(30) +
                                  f"{chi_w_ratio_error:.14e}".ljust(30) +
                                  f"{acceptance_rates_local_moves[0]:.14e}".ljust(50) +
                                  f"{acceptance_rates_local_moves[1]:.14e}".ljust(50) +
                                  f"{acceptance_rates_all_moves[0]:.14e}".ljust(50) +
                                  f"{acceptance_rates_all_moves[1]:.14e}".ljust(50) +
                                  f"{acceptance_rates_all_moves[2]:.14e}" + "\n")
            temperatures.append(temperature)
            chi_w_ratios.append(chi_w_ratio)

            temperature -= magnitude_of_temperature_increments

        print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
        pool.close()
        output_file.close()

    plt.plot(temperatures, chi_w_ratios, marker=".", markersize=5, color="k", linestyle='dashed')
    plt.xlabel(r"temperature, $1 / (\beta J)$", fontsize=15, labelpad=10)
    plt.ylabel(r"$\chi_{\rm w, local}$ / $\chi_{\rm w, all}$", fontsize=15, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
    plt.savefig(f"{output_directory}/prb_91_155412_fig_2_{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_"
                f"{algorithm_name.replace('-', '_')}.pdf", bbox_inches="tight")


if __name__ == "__main__":
    if len(sys.argv) > 2 or len(sys.argv) < 1:
        raise Exception("InterfaceError: provide at most one positional argument.  This is not required, but you may "
                        "provide the value False in order to choose the HXY version of the script (default value is "
                        "True, which chooses the elementary electrolyte).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be True / False to choose the elementary-electrolyte / HXY "
              "version of the script.")
        if not (sys.argv[1] == "True" or sys.argv[1] == "False"):
            raise Exception("InterfaceError: If provided, the single positional argument must be True or False.")
        main(eval(sys.argv[1]))
    else:
        print("Positional argument not provided.  The default value is True, which chooses the elementary-electrolyte "
              "version of the script.")
        main(True)
