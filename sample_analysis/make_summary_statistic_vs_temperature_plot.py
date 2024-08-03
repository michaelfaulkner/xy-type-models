from sample_getter import get_acceptance_rates, get_event_rate
from setup_scripts import check_initial_run_index, check_for_observable_vs_model_error, setup_pool
from markov_chain_diagnostics import get_sample_mean_and_error
import importlib
import matplotlib.pyplot as plt
import numpy as np
import os
import sample_getter
import sys
import time

# import additional modules; have to add the directory that contains run.py to sys.path
this_directory = os.path.dirname(os.path.abspath(__file__))
directory_containing_run_script = os.path.abspath(this_directory + "/../")
sys.path.insert(0, directory_containing_run_script)
run_script = importlib.import_module("run")


def main(config_file, observable_string):
    (algorithm_name, sample_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps, no_of_samples,
     temperatures, _, external_global_moves_string, no_of_runs, initial_run_index, max_no_of_cpus
     ) = run_script.get_config_data(config_file)
    check_for_observable_error(algorithm_name, observable_string)
    check_initial_run_index(initial_run_index)
    pool = setup_pool(no_of_runs, max_no_of_cpus)
    start_time = time.time()
    means, errors = get_means_and_errors(observable_string, algorithm_name, temperatures, no_of_sites,
                                         no_of_sites_string, no_of_equilibration_sweeps, no_of_samples,
                                         external_global_moves_string, sample_directory, sample_directory, no_of_runs,
                                         pool)
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")
    if no_of_runs > 1:
        pool.close()
    if not ((observable_string == "acceptance_rates") or (observable_string == "event_rate")):
        plt.errorbar(temperatures, means, errors, marker=".", markersize=5, color="k")
        plt.xlabel(r"temperature, $1 / (\beta J)$", fontsize=15, labelpad=10)
        plt.ylabel(f"{observable_string.replace('_', ' ')}", fontsize=15, labelpad=10)
        plt.tick_params(axis="both", which="major", labelsize=14, pad=10)
        plt.savefig(f"{sample_directory}/{observable_string}_vs_temperature_{algorithm_name.replace('-', '_')}_"
                    f"{external_global_moves_string}_{no_of_sites_string}.pdf", bbox_inches="tight")


def get_means_and_errors(observable_string, algorithm_name, temperatures, no_of_sites, no_of_sites_string,
                         no_of_equilibration_sweeps, no_of_samples, external_global_moves_string, output_directory,
                         sample_directory, no_of_runs, pool):
    output_file_string = (
        f"{output_directory}/{observable_string}_vs_temperature_{algorithm_name.replace('-', '_')}_"
        f"{external_global_moves_string}_{no_of_sites_string}_{no_of_runs}x{no_of_samples}_obs.tsv")
    try:
        with open(output_file_string, "r") as output_file:
            output_file_sans_header = np.array([np.fromstring(line, dtype=float, sep='\t') for line in output_file
                                                if not line.startswith('#')]).transpose()
            return output_file_sans_header[1], output_file_sans_header[2]
    except IOError:
        output_file = open(output_file_string, "w")
        if observable_string == "acceptance_rates":
            if no_of_runs == 1:
                acceptance_rates = get_acceptance_rates(sample_directory, 0)
            else:
                acceptance_rates = get_acceptance_rates(sample_directory + "/run_1", 0)
            if len(acceptance_rates) == 2:
                output_file.write("# temperature".ljust(30) + "final width of proposal interval".ljust(40) +
                                  "rotational acceptance rate" + "\n")
            elif len(acceptance_rates) == 3:
                if algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte":
                    output_file.write("# temperature".ljust(30) + "final width of proposal interval".ljust(40) +
                                      "acceptance rate (rotational moves)".ljust(40) +
                                      "acceptance rate (charge hops)" + "\n")
                else:
                    output_file.write("# temperature".ljust(30) + "final width of proposal interval".ljust(40) +
                                      "acceptance rate (rotational moves)".ljust(40) +
                                      "acceptance rate (external global moves)" + "\n")
            else:
                output_file.write("# temperature".ljust(30) + "final width of proposal interval".ljust(40) +
                                  "acceptance rate (field rotations)".ljust(40) +
                                  "acceptance rate (charge hops)".ljust(40) + "acceptance rate (global moves)" + "\n")
        elif observable_string == "event_rate":
            if no_of_runs == 1:
                event_rate = get_event_rate(sample_directory, 0)
            else:
                event_rate = get_event_rate(sample_directory + "/run_1", 0)
            if len(event_rate) == 1:
                output_file.write("# temperature".ljust(30) + "event rate (field rotations)" + "\n")
            else:
                output_file.write("# temperature".ljust(30) + "event rate (field rotations)".ljust(40) +
                                  "acceptance rate (external global moves)" + "\n")
        else:
            output_file.write("# temperature".ljust(30) + observable_string.ljust(35) + observable_string + " error" +
                              "\n")

        means, errors = [], []
        print(f"Processing {observable_string} of {algorithm_name}")
        for temperature_index, temperature in enumerate(temperatures):
            print(f"Temperature = {temperature:.4f}")
            get_sample_method = getattr(sample_getter, "get_" + observable_string)
            if observable_string == "acceptance_rates" or observable_string == "event_rate":
                if no_of_runs == 1:
                    acceptance_rates_or_event_rate = get_sample_method(sample_directory, temperature_index)
                else:
                    acceptance_rates_or_event_rate = np.mean([
                        get_sample_method(sample_directory + "/run_" + str(run_number), temperature_index)
                        for run_number in range(no_of_runs)], axis=0)
                if len(acceptance_rates_or_event_rate) == 1:
                    output_file.write(f"{temperature:.14e}".ljust(30) +
                                      f"{acceptance_rates_or_event_rate[0]:.14e}" + "\n")
                elif len(acceptance_rates_or_event_rate) == 2:
                    output_file.write(f"{temperature:.14e}".ljust(30) +
                                      f"{acceptance_rates_or_event_rate[0]:.14e}".ljust(40) +
                                      f"{acceptance_rates_or_event_rate[1]:.14e}" + "\n")
                elif len(acceptance_rates_or_event_rate) == 3:
                    output_file.write(f"{temperature:.14e}".ljust(30) +
                                      f"{acceptance_rates_or_event_rate[0]:.14e}".ljust(40) +
                                      f"{acceptance_rates_or_event_rate[1]:.14e}".ljust(40) +
                                      f"{acceptance_rates_or_event_rate[2]:.14e}" + "\n")
                else:
                    output_file.write(f"{temperature:.14e}".ljust(30) +
                                      f"{acceptance_rates_or_event_rate[0]:.14e}".ljust(40) +
                                      f"{acceptance_rates_or_event_rate[1]:.14e}".ljust(40) +
                                      f"{acceptance_rates_or_event_rate[2]:.14e}".ljust(40) +
                                      f"{acceptance_rates_or_event_rate[3]:.14e}" + "\n")
            else:
                if no_of_runs == 1:
                    sample_mean, sample_error = get_sample_mean_and_error(get_sample_method(
                        sample_directory, temperature, temperature_index, no_of_sites, no_of_equilibration_sweeps))
                else:
                    sample_means_and_errors = np.transpose(np.array(pool.starmap(get_sample_mean_and_error, [[
                        get_sample_method(f"{sample_directory}/run_{run_number}", temperature, temperature_index,
                                          no_of_sites, no_of_equilibration_sweeps)]
                        for run_number in range(no_of_runs)])))
                    """comment out the above (replacing w/below) if having problems w/pool.starmap() on cluster"""
                    """
                    sample_means_and_errors = np.transpose(np.array([get_sample_mean_and_error(get_sample_method(
                        f"{sample_directory}/run_{run_number}", temperature, temperature_index, no_of_sites,
                        no_of_equilibration_sweeps)) for run_number in range(no_of_runs)]))
                    """
                    sample_mean = np.mean(sample_means_and_errors[0])
                    sample_error = np.linalg.norm(sample_means_and_errors[1])
                output_file.write(f"{temperature:.14e}".ljust(30) + f"{sample_mean:.14e}".ljust(35) +
                                  f"{sample_error:.14e}" + "\n")
                means.append(sample_mean), errors.append(sample_error)
        output_file.close()
        return means, errors


def check_for_observable_error(algorithm_name, observable_string):
    # Check to ensure that either acceptance_rates, event_rate or a one-dimensional observable has been given as the
    # second positional argument
    if (observable_string != "acceptance_rates" and observable_string != "event_rate"
            and observable_string != "potential" and observable_string != "specific_heat" 
            and observable_string != "magnetisation_norm" and observable_string != "magnetisation_phase"
            and observable_string != "rotated_magnetisation_phase" and observable_string != "magnetisation_squared"
            and observable_string != "magnetic_susceptibility" and observable_string != "relative_magnetisation_norm"
            and observable_string != "inverse_vacuum_permittivity" and observable_string != "helicity_modulus"
            and observable_string != "hxy_internal_twist_susceptibility"
            and observable_string != "xy_twist_relaxation_susceptibility"
            and observable_string != "xy_emergent_field_zero_mode_susceptibility"
            and observable_string != "xy_global_defect_susceptibility"
            and observable_string != "potential_minimising_twist_susceptibility"
            and observable_string != "inverse_permittivity" and observable_string != "topological_susceptibility"):
        print("ConfigurationError: Give one of acceptance_rates, event_rate, potential, specific_heat, "
              "magnetisation_norm, magnetisation_phase, rotated_magnetisation_phase, magnetisation_squared, "
              "magnetic_susceptibility, relative_magnetisation_norm, inverse_vacuum_permittivity, helicity_modulus, "
              "hxy_internal_twist_susceptibility, xy_twist_relaxation_susceptibility, "
              "xy_emergent_field_zero_mode_susceptibility, xy_global_defect_susceptibility, "
              "potential_minimising_twist_susceptibility, inverse_permittivity or topological_susceptibility as the "
              "second positional argument.")
        raise SystemExit
    # Raise an error if acceptance_rates has been given as the second positional argument for sample generated by an 
    # ECMC algorithm
    if (algorithm_name == "xy-ecmc" or algorithm_name == "hxy-ecmc") and observable_string == "acceptance_rates":
        print("ConfigurationError: These samples were generated by an event-chain Monte Carlo algorithm: rather than "
              "acceptance_rates, give event_rate as the second positional argument.")
        raise SystemExit
    # Raise an error if event_rate has been given as the second positional argument for sample generated by a
    # Metropolis algorithm
    if ((algorithm_name == "xy-uniform-noise-metropolis" or algorithm_name == "hxy-uniform-noise-metropolis" or
         algorithm_name == "xy-gaussian-noise-metropolis" or algorithm_name == "hxy-gaussian-noise-metropolis" or
         algorithm_name == "elementary-electrolyte" or algorithm_name == "multivalued-electrolyte") and
            observable_string == "event_rate"):
        print("ConfigurationError: These samples were generated by a Metropolis algorithm: rather than "
              "event_rate, give acceptance_rates as the second positional argument.")
        raise SystemExit
    check_for_observable_vs_model_error(algorithm_name, observable_string)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose summary statistics you wish to calculate in the second and "
                        "third positions, respectively.")
    main(sys.argv[1], sys.argv[2])
