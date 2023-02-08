from markov_chain_diagnostics import get_sample_mean_and_error
from sample_getter import get_acceptance_rates
import importlib
import math
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import sample_getter
import sys
import time

# import run script - have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def main(model):
    if model == "electrolyte":
        approx_transition_temperature = 1.351
        config_file_local_moves = "config_files/prb_91_155412_fig_2/local_moves.txt"
        config_file_all_moves = "config_files/prb_91_155412_fig_2/all_moves.txt"
        get_sample_methods = [getattr(sample_getter, "get_topological_susceptibility")]
    elif model == "hxy":
        approx_transition_temperature = 1.351
        config_file_local_moves = "config_files/prb_91_155412_fig_2_hxy/local_moves.txt"
        config_file_all_moves = "config_files/prb_91_155412_fig_2_hxy/all_moves.txt"
        get_sample_methods = [getattr(sample_getter, "get_potential_minimising_twist_susceptibility"),
                              getattr(sample_getter, "get_helicity_modulus"),
                              getattr(sample_getter, "get_hxy_internal_twist_susceptibility")]
    elif model == "xy":
        approx_transition_temperature = 0.887
        config_file_local_moves = "config_files/prb_91_155412_fig_2_xy/local_moves.txt"
        config_file_all_moves = "config_files/prb_91_155412_fig_2_xy/all_moves.txt"
        get_sample_methods = [getattr(sample_getter, "get_potential_minimising_twist_susceptibility"),
                              getattr(sample_getter, "get_helicity_modulus")]
    else:
        raise Exception("InterfaceError: If provided, the single positional argument must be electrolyte, hxy or xy.")
    (algorithm_name, output_directory_local_moves, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps, _,
     temperatures, _, _, no_of_jobs, _, max_no_of_cpus) = run_script.get_config_data(config_file_local_moves)
    output_directory_all_moves = run_script.get_config_data(config_file_all_moves)[1]

    try:
        get_acceptance_rates(f"{output_directory_local_moves}/job_0", 0)
    except OSError:
        raise Exception(f"Local-move simulations have not been run - enter 'python run.py {config_file_local_moves}' "
                        f"in the top directory.")
    try:
        get_acceptance_rates(f"{output_directory_all_moves}/job_0", 0)
    except OSError:
        raise Exception(f"Local-move simulations have not been run - enter 'python run.py {config_file_all_moves}' in "
                        f"the top directory.")

    output_directory = output_directory_local_moves.replace("/local_moves", "")
    reduced_temperatures = [temperature / approx_transition_temperature for temperature in temperatures]

    accept_rates_file_string = (f"{output_directory}/acceptance_rates_{algorithm_name.replace('-', '_')}_"
                                f"{no_of_sites_string}.tsv")
    output_file_strings = [f"{output_directory}/{algorithm_name.replace('-', '_')}_{no_of_sites_string}_topological_"
                           f"susceptibility_ratio.tsv"]
    figure_file_strings = [f"{output_directory}/{algorithm_name.replace('-', '_')}_{no_of_sites_string}_topological_"
                           f"susceptibility_ratio.pdf"]
    if model != "electrolyte":
        output_file_strings.append(f"{output_directory}/{algorithm_name.replace('-', '_')}_{no_of_sites_string}_"
                                   f"helicity_ratio.tsv")
        figure_file_strings.append(f"{output_directory}/{algorithm_name.replace('-', '_')}_{no_of_sites_string}_"
                                   f"helicity_ratio.pdf")
    if model == "hxy":
        output_file_strings.append(f"{output_directory}/{algorithm_name.replace('-', '_')}_{no_of_sites_string}_"
                                   f"internal_twist_susceptibility_ratio.tsv")
        figure_file_strings.append(f"{output_directory}/{algorithm_name.replace('-', '_')}_{no_of_sites_string}_"
                                   f"internal_twist_susceptibility_ratio.pdf")

    no_of_cpus = mp.cpu_count()
    pool = mp.Pool(no_of_cpus)
    start_time = time.time()

    accept_rates_file = open(accept_rates_file_string, "w")
    if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
        accept_rates_file.write("# temperature".ljust(30) +
                                "final width of proposal interval (local only)".ljust(50) +
                                "rotational acceptance rate (local only)".ljust(50) +
                                "charge-hop acceptance rate (local only)".ljust(50) +
                                "final width of proposal interval (all moves)".ljust(50) +
                                "rotational acceptance rate (all moves)".ljust(50) +
                                "charge-hop acceptance rate (all moves)".ljust(50) +
                                "acceptance rate (external global moves)" + "\n")
    else:
        accept_rates_file.write("# temperature".ljust(30) +
                                "final width of proposal interval (local only)".ljust(50) +
                                "rotational acceptance rate (local only)".ljust(50) +
                                "final width of proposal interval (all moves)".ljust(50) +
                                "rotational acceptance rate (all moves)".ljust(50) +
                                "acceptance rate (external global moves)" + "\n")
    for temperature_index, temperature in enumerate(temperatures):
        acceptance_rates_local_moves = np.mean(
            [get_acceptance_rates(output_directory_local_moves + f"/job_{job_number}", temperature_index) for
             job_number in range(no_of_jobs)], axis=0)
        acceptance_rates_all_moves = np.mean(
            [get_acceptance_rates(output_directory_all_moves + f"/job_{job_number}", temperature_index) for
             job_number in range(no_of_jobs)], axis=0)
        if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
            accept_rates_file.write(f"{temperature:.14e}".ljust(30) +
                                    f"{acceptance_rates_local_moves[0]:.14e}".ljust(50) +
                                    f"{acceptance_rates_local_moves[1]:.14e}".ljust(50) +
                                    f"{acceptance_rates_local_moves[2]:.14e}".ljust(50) +
                                    f"{acceptance_rates_all_moves[0]:.14e}".ljust(50) +
                                    f"{acceptance_rates_all_moves[1]:.14e}".ljust(50) +
                                    f"{acceptance_rates_all_moves[2]:.14e}".ljust(50) +
                                    f"{acceptance_rates_all_moves[3]:.14e}" + "\n")
        else:
            accept_rates_file.write(f"{temperature:.14e}".ljust(30) +
                                    f"{acceptance_rates_local_moves[0]:.14e}".ljust(50) +
                                    f"{acceptance_rates_local_moves[1]:.14e}".ljust(50) +
                                    f"{acceptance_rates_all_moves[0]:.14e}".ljust(50) +
                                    f"{acceptance_rates_all_moves[1]:.14e}".ljust(50) +
                                    f"{acceptance_rates_all_moves[2]:.14e}" + "\n")
    accept_rates_file.close()

    for observable_index, output_file_string in enumerate(output_file_strings):
        # in the following code, chi represents the observable corresponding to observable_index
        get_sample_method = get_sample_methods[observable_index]
        try:
            with open(output_file_string, "r") as output_file:
                output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                    if not line.startswith('#')]).transpose()
                chi_ratios = output_file_sans_header[1]
        except IOError:
            output_file = open(output_file_string, "w")
            output_file.write("# *** NB, chi represents the observable in the file name ***" + "\n")
            output_file.write("# temperature".ljust(30) + "chi ratio".ljust(30) + "chi ratio error".ljust(30) +
                              "chi (local only)".ljust(30) + "chi error (local only)".ljust(30) +
                              "chi (all moves)".ljust(30) + "chi error (all moves)".ljust(30) + "\n")
            chi_ratios = []
            for temperature_index, temperature in enumerate(temperatures):
                print(f"Temperature = {temperature:.4f}")
                sample_means_and_errors_local_moves = np.transpose(np.array(pool.starmap(get_sample_mean_and_error, [[
                    get_sample_method(
                        output_directory_local_moves + f"/job_{job_number}", temperature, temperature_index,
                        no_of_sites, no_of_equilibration_sweeps)] for job_number in range(no_of_jobs)])))
                sample_means_and_errors_all_moves = np.transpose(np.array(pool.starmap(get_sample_mean_and_error, [[
                    get_sample_method(output_directory_all_moves + f"/job_{job_number}", temperature, temperature_index,
                                      no_of_sites, no_of_equilibration_sweeps)] for job_number in range(no_of_jobs)])))

                sample_mean_local_moves = np.mean(sample_means_and_errors_local_moves[0])
                sample_error_local_moves = np.linalg.norm(sample_means_and_errors_local_moves[1])
                sample_mean_all_moves = np.mean(sample_means_and_errors_all_moves[0])
                sample_error_all_moves = np.linalg.norm(sample_means_and_errors_all_moves[1])
                chi_ratio = sample_mean_local_moves / sample_mean_all_moves
                chi_ratio_error = math.sqrt((sample_error_local_moves / sample_mean_all_moves) ** 2 +
                                            (sample_mean_local_moves * sample_error_all_moves /
                                             sample_mean_all_moves ** 2) ** 2)

                output_file.write(f"{temperature:.14e}".ljust(30) + f"{chi_ratio:.14e}".ljust(30) +
                                  f"{chi_ratio_error:.14e}".ljust(30) + f"{sample_mean_local_moves:.14e}".ljust(30) +
                                  f"{sample_error_local_moves:.14e}".ljust(30) +
                                  f"{sample_mean_all_moves:.14e}".ljust(30) +
                                  f"{sample_error_all_moves:.14e}".ljust(30) + "\n")
                chi_ratios.append(chi_ratio)
            output_file.close()

        plt.plot(reduced_temperatures, chi_ratios, marker=".", markersize=5, color="k", linestyle='dashed')
        plt.xlabel(r"reduced temperature, $\beta_{\rm c} / \beta$", fontsize=15, labelpad=10)
        plt.ylabel(r"$\chi_{\rm{local}}$ / $\chi_{\rm{all}}$", fontsize=15, labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
        plt.ylim([-0.05, 1.05])
        plt.savefig(figure_file_strings[observable_index], bbox_inches="tight")
        plt.clf()

    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    pool.close()


if __name__ == "__main__":
    if len(sys.argv) > 2 or len(sys.argv) < 1:
        raise Exception("InterfaceError: provide at most one positional argument.  You may provide the value "
                        "electrolyte / hxy / xy in order to choose the elementary-electrolyte / HXY / XY version of the"
                        " script (default value is electrolyte).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be electrolyte / hxy / xy to choose the "
              "elementary-electrolyte / HXY / XY version of the script.")
        main(str(sys.argv[1]))
    else:
        print("Positional argument not provided.  The default value is electrolyte, which chooses the "
              "elementary-electrolyte version of the script.")
        main("electrolyte")
