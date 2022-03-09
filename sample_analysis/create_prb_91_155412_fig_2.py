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


def main(model):
    if model == "electrolyte":
        config_file_local_moves = "config_files/prb_91_155412_fig_2/local_moves.txt"
        config_file_all_moves = "config_files/prb_91_155412_fig_2/all_moves.txt"
        get_sample_method = getattr(sample_getter, "get_topological_susceptibility")
    elif model == "hxy":
        config_file_local_moves = "config_files/prb_91_155412_fig_2_hxy/local_moves.txt"
        config_file_all_moves = "config_files/prb_91_155412_fig_2_hxy/all_moves.txt"
        get_sample_method = getattr(sample_getter, "get_potential_minimising_twist_susceptibility")
    elif model == "hxy_literal":
        config_file_local_moves = "config_files/prb_91_155412_fig_2_hxy_literal/local_moves.txt"
        config_file_all_moves = "config_files/prb_91_155412_fig_2_hxy_literal/all_moves.txt"
        get_sample_method = getattr(sample_getter, "get_hxy_topological_susceptibility")
    elif model == "xy":
        config_file_local_moves = "config_files/prb_91_155412_fig_2_xy/local_moves.txt"
        config_file_all_moves = "config_files/prb_91_155412_fig_2_xy/all_moves.txt"
        get_sample_method = getattr(sample_getter, "get_potential_minimising_twist_susceptibility")
    else:
        raise Exception("InterfaceError: If provided, the single positional argument must be electrolyte, hxy, "
                        "hxy_literal or xy.")
    (algorithm_name, output_directory_local_moves, no_of_sites, no_of_equilibration_sweeps, _, temperatures, _, _,
     no_of_jobs, max_no_of_cpus) = run_script.get_config_data(config_file_local_moves)
    output_directory_all_moves = run_script.get_config_data(config_file_all_moves)[1]
    try:
        sample_getter.get_acceptance_rates(f"{output_directory_local_moves}/job_0", 0)
    except OSError:
        raise Exception(f"Local-move simulations have not been run - enter 'python run.py {config_file_local_moves}' "
                        f"in the top directory.")
    try:
        sample_getter.get_acceptance_rates(f"{output_directory_all_moves}/job_0", 0)
    except OSError:
        raise Exception(f"Local-move simulations have not been run - enter 'python run.py {config_file_all_moves}' in "
                        f"the top directory.")
    output_directory = output_directory_local_moves.replace("/local_moves", "")

    if model == "hxy_literal":
        output_file_string = (f"{output_directory}/prb_91_155412_fig_2_literal_{algorithm_name.replace('-', '_')}_"
                              f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites.tsv")
        accept_rates_file_string = (f"{output_directory}/acceptance_rates_literal_{algorithm_name.replace('-', '_')}_"
                                    f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites.tsv")
        figure_file_string = (f"{output_directory}/prb_91_155412_fig_2_literal_{algorithm_name.replace('-', '_')}_"
                              f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites.pdf")
    else:
        output_file_string = (f"{output_directory}/prb_91_155412_fig_2_{algorithm_name.replace('-', '_')}_"
                              f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites.tsv")
        accept_rates_file_string = (f"{output_directory}/acceptance_rates_{algorithm_name.replace('-', '_')}_"
                                    f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites.tsv")
        figure_file_string = (f"{output_directory}/prb_91_155412_fig_2_{algorithm_name.replace('-', '_')}_"
                              f"{int(no_of_sites ** 0.5)}x{int(no_of_sites ** 0.5)}_sites.pdf")

    try:
        with open(output_file_string, "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            chi_ratios = output_file_sans_header[1]
    except (IOError, IndexError) as _:
        output_file = open(output_file_string, "w")
        accept_rates_file = open(accept_rates_file_string, "w")
        output_file.write("# temperature".ljust(30) + "chi ratio".ljust(30) + "chi ratio error".ljust(30) +
                          "chi (local only)".ljust(30) + "chi error (local only)".ljust(30) +
                          "chi (all moves)".ljust(30) + "chi error (all moves)".ljust(30) + "\n")
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

        no_of_cpus = mp.cpu_count()
        pool = mp.Pool(no_of_cpus)
        chi_ratios = []
        start_time = time.time()

        for temperature_index, temperature in enumerate(temperatures):
            print(f"Temperature = {temperature:.4f}")
            acceptance_rates_local_moves = np.mean(
                [sample_getter.get_acceptance_rates(output_directory_local_moves + f"/job_{job_number}",
                                                    temperature_index) for job_number in range(no_of_jobs)], axis=0)
            acceptance_rates_all_moves = np.mean(
                [sample_getter.get_acceptance_rates(output_directory_all_moves + f"/job_{job_number}",
                                                    temperature_index) for job_number in range(no_of_jobs)], axis=0)

            sample_means_and_errors_local_moves = np.transpose(
                np.array(pool.starmap(markov_chain_diagnostics.get_sample_mean_and_error, [[get_sample_method(
                        output_directory_local_moves + f"/job_{job_number}", temperature, temperature_index,
                        no_of_sites)[no_of_equilibration_sweeps:]] for job_number in range(no_of_jobs)])))
            sample_means_and_errors_all_moves = np.transpose(
                np.array(pool.starmap(markov_chain_diagnostics.get_sample_mean_and_error, [[get_sample_method(
                        output_directory_all_moves + f"/job_{job_number}", temperature, temperature_index,
                        no_of_sites)[no_of_equilibration_sweeps:]] for job_number in range(no_of_jobs)])))

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
            chi_ratios.append(chi_ratio)
        print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
        pool.close()
        output_file.close()

    plt.plot(temperatures, chi_ratios, marker=".", markersize=5, color="k", linestyle='dashed')
    plt.xlabel(r"temperature, $1 / (\beta J)$", fontsize=15, labelpad=10)
    plt.ylabel(r"$\chi_{\rm{local}}$ / $\chi_{\rm{all}}$", fontsize=15, labelpad=10)
    plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
    plt.savefig(figure_file_string, bbox_inches="tight")


if __name__ == "__main__":
    if len(sys.argv) > 2 or len(sys.argv) < 1:
        raise Exception("InterfaceError: provide at most one positional argument.  You may provide the value "
                        "electrolyte / hxy / hxy_literal / xy in order to choose the elementary-electrolyte / HXY / "
                        "HXY with topological susceptibility / XY version of the script (default value is "
                        "electrolyte).")
    if len(sys.argv) == 2:
        print("One positional argument provided.  This must be electrolyte / hxy / hxy_literal / xy to choose the "
              "elementary-electrolyte / HXY / HXY with topological susceptibility / XY version of the script.")
        main(str(sys.argv[1]))
    else:
        print("Positional argument not provided.  The default value is electrolyte, which chooses the "
              "elementary-electrolyte version of the script.")
        main("electrolyte")
