import matplotlib.pyplot as plt
import math
import numpy as np
import os


def main():
    data_directory = "permanent_data/spin_config_plots"
    output_directory = "spin_config_plots"
    xy_spin_phases = np.load(f"{data_directory}/xy_snapshot_temp_eq_1_spin_phases.npy")
    xy_spin_phases_zero_temp = np.load(f"{data_directory}/xy_snapshot_temp_eq_1_spin_phases_zero_temp.npy")
    xy_spin_phases_hxy_zero_temp = np.load(f"{data_directory}/xy_snapshot_temp_eq_1_spin_phases_hxy_zero_temp.npy")
    hxy_spin_phases = np.load(f"{data_directory}/hxy_snapshot_temp_eq_1.5_spin_phases.npy")
    hxy_spin_phases_zero_temp = np.load(f"{data_directory}/hxy_snapshot_temp_eq_1.5_spin_phases_zero_temp.npy")

    make_single_plot(xy_spin_phases, xy_spin_phases_zero_temp, output_directory,
                     "decomposed_spin_configs_xy_temp_eq_1.pdf")
    make_single_plot(xy_spin_phases, xy_spin_phases_hxy_zero_temp, output_directory,
                     "decomposed_spin_configs_xy_with_hxy_minimisation_temp_eq_1.pdf")
    make_single_plot(hxy_spin_phases, hxy_spin_phases_zero_temp, output_directory,
                     "decomposed_spin_configs_hxy_temp_eq_1.5.pdf")


def make_single_plot(spin_phases, zero_temp_spin_phases, output_directory, output_file_string):
    lattice_coordinates = [[index // 20 + 1 for index in range(len(spin_phases))],
                           [index % 20 + 1 for index in range(len(spin_phases))]]
    untwisted_spin_field = [
        [np.cos(spin_phase) * math.cos(math.pi * (index % 20) / 10) +
         np.sin(spin_phase) * math.sin(math.pi * (index % 20) / 10) for index, spin_phase in
         enumerate(zero_temp_spin_phases)],
        [np.sin(spin_phase) * math.cos(math.pi * (index % 20) / 10) -
         np.cos(spin_phase) * math.sin(math.pi * (index % 20) / 10) for index, spin_phase in
         enumerate(zero_temp_spin_phases)]]
    global_twist_field = [[math.sin(math.pi * (index % 20) / 10) for index in range(len(spin_phases))],
                          [- math.cos(math.pi * (index % 20) / 10) for index in range(len(spin_phases))]]

    fig, axes = plt.subplots(1, 4, figsize=(16.0, 4.0))
    fig.tight_layout(w_pad=-12)
    axes[0].quiver(lattice_coordinates[0], lattice_coordinates[1], np.cos(spin_phases), np.sin(spin_phases), scale=25,
                   width=0.005, color="r", pivot="mid")
    axes[1].quiver(lattice_coordinates[0], lattice_coordinates[1], np.cos(zero_temp_spin_phases),
                   np.sin(zero_temp_spin_phases), scale=25, width=0.005, color="r", pivot="mid")
    axes[2].quiver(lattice_coordinates[0], lattice_coordinates[1], untwisted_spin_field[0], untwisted_spin_field[1],
                   scale=25, width=0.005, color="r", pivot="mid")
    axes[3].quiver(lattice_coordinates[0], lattice_coordinates[1], global_twist_field[0], global_twist_field[1],
                   scale=25,
                   width=0.005, color="r", pivot="mid")

    [axis.axis("square") for axis in axes]
    [axis.xaxis.set_tick_params(labelbottom=False) for axis in axes]
    [axis.yaxis.set_tick_params(labelleft=False) for axis in axes]
    [axis.set_xticks([]) for axis in axes], [axis.set_yticks([]) for axis in axes]
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes]
    letter_heights = 0.975
    fig.text(0.07125, letter_heights, "a", fontsize=25, weight='bold')
    fig.text(0.289, letter_heights, "b", fontsize=25, weight='bold')
    fig.text(0.509, letter_heights, "c", fontsize=25, weight='bold')
    fig.text(0.7285, letter_heights, "d", fontsize=25, weight='bold')
    os.system(f"mkdir -p {output_directory}")
    fig.savefig(f"{output_directory}/{output_file_string}", bbox_inches="tight")
    fig.clear()


if __name__ == "__main__":
    main()
