from polyspectra import get_power_trispectrum, get_power_trispectrum_as_defined
from setup_scripts import reverse_enumerate, setup_polyspectra_script
import matplotlib.pyplot as plt
import sys
import time


def main(config_file, observable_string, no_of_trispectrum_octaves=3, trispectrum_base_period_shift=1):
    (algorithm_name, output_directory, no_of_sites, no_of_sites_string, no_of_equilibration_sweeps, temperatures,
     external_global_moves_string, no_of_runs, pool) = setup_polyspectra_script(config_file, observable_string)

    start_time = time.time()
    for temperature_index, temperature in reverse_enumerate(temperatures):
        print(f"Temperature = {temperature:.4f}")

        power_trispectrum = get_power_trispectrum(
            algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
            no_of_sites_string, external_global_moves_string, no_of_runs, pool, no_of_trispectrum_octaves,
            trispectrum_base_period_shift, no_of_equilibration_sweeps)
        power_trispectrum_as_defined = get_power_trispectrum_as_defined(
            algorithm_name, observable_string, output_directory, temperature, temperature_index, no_of_sites,
            no_of_sites_string, external_global_moves_string, no_of_runs, pool, no_of_trispectrum_octaves,
            trispectrum_base_period_shift, no_of_equilibration_sweeps)

        # normalise each estimator of power trispectrum with respect to their low-frequency values
        power_trispectrum[2] = [spectrum / spectrum[0] for spectrum in power_trispectrum[2]]
        power_trispectrum_as_defined[2] = [spectrum / spectrum[0] for spectrum in power_trispectrum_as_defined[2]]

        figure, axis = plt.subplots(2, 2, figsize=(10, 10))
        [axis[1, index].set_xlabel(r"frequency, $f$ $(t^{-1})$", fontsize=10, labelpad=10) for index in range(2)]
        plt.tick_params(axis="both", which="major", labelsize=10, pad=10)

        for index in range(2):
            if index == 0:
                axis[index, 0].loglog(power_trispectrum_as_defined[1], power_trispectrum_as_defined[2][index],
                                      color='blue', label="estimator as defined")
                axis[index, 0].loglog(power_trispectrum[1], power_trispectrum[2][index], color='red',
                                      label="shortcut estimator")

                axis[index, 1].loglog(power_trispectrum[1], power_trispectrum[2][index], color='red',
                                      label="shortcut estimator")
                axis[index, 1].loglog(power_trispectrum_as_defined[1], power_trispectrum_as_defined[2][index],
                                      color='blue', label="estimator as defined")

                axis[index, 0].set_ylabel(fr"$|S_X^3 \left( f, f_0' \right)|$ / $|S_X^3 \left( f_0, f_0' \right)|$, "
                                          fr"$f_0' = {power_trispectrum[0][0]:.2e}$", fontsize=10, labelpad=10)
            else:
                axis[index, 0].loglog(power_trispectrum_as_defined[1],
                                      power_trispectrum_as_defined[2][len(power_trispectrum_as_defined[2]) - 1],
                                      color='blue', label="estimator as defined")
                axis[index, 0].loglog(power_trispectrum[1],
                                      power_trispectrum[2][len(power_trispectrum[2]) - 1], color='red',
                                      label="shortcut estimator")

                axis[index, 1].loglog(power_trispectrum[1],
                                      power_trispectrum[2][len(power_trispectrum[2]) - 1], color='red',
                                      label="shortcut estimator")
                axis[index, 1].loglog(power_trispectrum_as_defined[1],
                                      power_trispectrum_as_defined[2][len(power_trispectrum_as_defined[2]) - 1],
                                      color='blue', label="estimator as defined")

                axis[index, 0].set_ylabel(fr"$|S_X^3 \left( f, {2 ** no_of_trispectrum_octaves} f_0' \right)|$ / "
                                          fr"$|S_X^3 \left(f_0, {2 ** no_of_trispectrum_octaves} f_0' \right)|$, "
                                          fr"$f_0' = {power_trispectrum[0][0]:.2e}$", fontsize=10, labelpad=10)

        figure.tight_layout()
        trispectrum_legend = (axis[0, 0].legend(loc="lower left", fontsize=10),
                              axis[0, 1].legend(loc="lower left", fontsize=10))
        for legend in trispectrum_legend:
            legend.get_frame().set_edgecolor("k")
            legend.get_frame().set_lw(1.5)
        figure.savefig(f"{output_directory}/{observable_string}_compare_trispectrum_estimators_"
                       f"{algorithm_name.replace('-', '_')}_{external_global_moves_string}_{no_of_sites_string}_temp_eq"
                       f"_{temperature:.4f}_max_shift_eq_{2 ** no_of_trispectrum_octaves}_x_"
                       f"{trispectrum_base_period_shift}_delta_t.pdf", bbox_inches="tight")
        figure.clf()
    print(f"Sample analysis complete.  Total runtime = {time.time() - start_time:.2e} seconds.")

    if no_of_runs > 1:
        pool.close()


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        raise Exception("InterfaceError: Two positional arguments required - give the configuration-file location and "
                        "the string of the observable whose power trispectrum you wish to estimate (in the first and "
                        "second positions, respectively).  In addition, you may provide no_of_trispectrum_octaves "
                        "(default value is 3) and trispectrum_base_period_shift (default value is 1) in the third and "
                        "fourth positions, respectively.")
    if len(sys.argv) == 3:
        print("Two positional arguments provided.  The first / second must be the location of the configuration file / "
              "the string of the observable whose power trispectrum you wish to estimate.  In addition, you may provide"
              " no_of_trispectrum_octaves (default value is 3) and trispectrum_base_period_shift (default value is 1) "
              "in the third and fourth positions, respectively.")
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        print("Three positional arguments provided.  The first / second / third must be the location of the "
              "configuration file / the string of the observable whose power trispectrum you wish to estimate / "
              "no_of_trispectrum_octaves.  In addition, you may provide trispectrum_base_period_shift (default value "
              "is 1) in the fourth position.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    elif len(sys.argv) == 5:
        print("Four positional arguments provided.  The first / second / third / fourth must be the location of the "
              "configuration file / the string of the observable whose power trispectrum you wish to estimate / "
              "no_of_trispectrum_octaves / trispectrum_base_period_shift.")
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
