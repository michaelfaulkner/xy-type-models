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
    hxy_spin_phases = np.load(f"{data_directory}/hxy_snapshot_temp_eq_1_point_5_spin_phases.npy")
    hxy_spin_phases_zero_temp = np.load(f"{data_directory}/hxy_snapshot_temp_eq_1_point_5_spin_phases_zero_temp.npy")
    make_set_of_plots(xy_spin_phases, xy_spin_phases_zero_temp, output_directory, "xy_temp_eq_1")
    make_set_of_plots(xy_spin_phases, xy_spin_phases_hxy_zero_temp, output_directory,
                      "xy_with_hxy_minimisation_temp_eq_1")
    make_set_of_plots(hxy_spin_phases, hxy_spin_phases_zero_temp, output_directory, "hxy_temp_eq_1_point_5")


def make_set_of_plots(spin_phases, zero_temp_spin_phases, output_directory, output_file_string):
    lattice_length = int(len(spin_phases) ** 0.5)
    lattice_coordinates = np.array([[index % lattice_length + 1 for index in range(len(spin_phases))],
                                    [index // lattice_length + 1 for index in range(len(spin_phases))]])

    spin_wave_spin_phases = spin_phases - zero_temp_spin_phases - 0.5 * math.pi  # pi / 2 rotation improves presentation
    untwisted_spin_phases = [spin_phase - 2.0 * math.pi * (index // lattice_length) / lattice_length
                             for index, spin_phase in enumerate(zero_temp_spin_phases)]
    global_twist_spin_phases = [2.0 * math.pi * (index // lattice_length) / lattice_length
                                for index in range(len(zero_temp_spin_phases))]

    untwisted_vector_spin_field = [
        [np.cos(spin_phase) * math.cos(2.0 * math.pi * (index // lattice_length) / lattice_length) +
         np.sin(spin_phase) * math.sin(2.0 * math.pi * (index // lattice_length) / lattice_length)
         for index, spin_phase in enumerate(zero_temp_spin_phases)],
        [np.sin(spin_phase) * math.cos(2.0 * math.pi * (index // lattice_length) / lattice_length) -
         np.cos(spin_phase) * math.sin(2.0 * math.pi * (index // lattice_length) / lattice_length)
         for index, spin_phase in enumerate(zero_temp_spin_phases)]]
    global_twist_vector_spin_field = [
        [math.sin(2.0 * math.pi * (index // lattice_length) / lattice_length) for index in range(len(spin_phases))],
        [- math.cos(2.0 * math.pi * (index // lattice_length) / lattice_length) for index in range(len(spin_phases))]]

    zeros = [0.0] * (lattice_length ** 2)
    electric_field_total = get_emergent_field_from_spin_phases(spin_phases, lattice_length)
    electric_field_irrotational = get_emergent_field_from_spin_phases(zero_temp_spin_phases, lattice_length)
    electric_field_poisson = get_emergent_field_from_spin_phases(untwisted_spin_phases, lattice_length)
    electric_field_top_sector = get_emergent_field_from_spin_phases(global_twist_spin_phases, lattice_length)
    electric_field_auxiliary = get_emergent_field_from_spin_phases(spin_wave_spin_phases, lattice_length)

    fig, axes = plt.subplots(1, 4, figsize=(16.0, 4.0))
    fig.tight_layout(w_pad=-12)
    plot_total_and_coulomb_spin_fields_of_single_row_fig(axes, lattice_coordinates, spin_phases, zero_temp_spin_phases)
    axes[2].quiver(lattice_coordinates[0], lattice_coordinates[1], untwisted_vector_spin_field[0],
                   untwisted_vector_spin_field[1], scale=25, width=0.005, color="r", pivot="mid")
    axes[3].quiver(lattice_coordinates[0], lattice_coordinates[1], global_twist_vector_spin_field[0],
                   global_twist_vector_spin_field[1], scale=25, width=0.005, color="r", pivot="mid")

    setup_axes_of_single_row_fig(axes)
    letter_heights = 0.975
    fig.text(0.07125, letter_heights, "(a)", fontsize=25)
    fig.text(0.289, letter_heights, "(b)", fontsize=25)
    fig.text(0.509, letter_heights, "(c)", fontsize=25)
    fig.text(0.7285, letter_heights, "(d)", fontsize=25)
    os.system(f"mkdir -p {output_directory}")
    fig.savefig(f"{output_directory}/coulomb_poisson_global_twist_decomposition_{output_file_string}.pdf",
                bbox_inches="tight")
    fig.clear()

    fig, axes = plt.subplots(2, 4, figsize=(16.5, 8.0))
    fig.tight_layout(w_pad=-12)
    plot_total_and_coulomb_spin_fields_of_double_row_fig(axes, lattice_coordinates, spin_phases, zero_temp_spin_phases)
    axes[0, 2].quiver(lattice_coordinates[0], lattice_coordinates[1], untwisted_vector_spin_field[0],
                      untwisted_vector_spin_field[1], scale=25, width=0.005, color="r", pivot="mid")
    axes[0, 3].quiver(lattice_coordinates[0], lattice_coordinates[1], global_twist_vector_spin_field[0],
                      global_twist_vector_spin_field[1], scale=25, width=0.005, color="r", pivot="mid")

    plot_total_and_coulomb_emergent_fields(axes, lattice_coordinates, electric_field_total, electric_field_irrotational,
                                           zeros)
    axes[1, 2].quiver(lattice_coordinates[0], lattice_coordinates[1] - 0.5, electric_field_poisson[0], zeros,
                      scale=40, width=0.002, color="k", pivot="mid")
    axes[1, 2].quiver(lattice_coordinates[0] - 0.5, lattice_coordinates[1], zeros, electric_field_poisson[1],
                      scale=40, width=0.002, color="k", pivot="mid")
    axes[1, 3].quiver(lattice_coordinates[0], lattice_coordinates[1] - 0.5, electric_field_top_sector[0], zeros,
                      scale=10, width=0.002, color="k", pivot="mid")
    axes[1, 3].quiver(lattice_coordinates[0] - 0.5, lattice_coordinates[1], zeros, electric_field_top_sector[1],
                      scale=10, width=0.002, color="k", pivot="mid")

    setup_axes_of_double_row_fig(axes)
    letter_heights = 0.99
    fig.text(0.07125, letter_heights, "(a)", fontsize=25)
    fig.text(0.289, letter_heights, "(b)", fontsize=25)
    fig.text(0.509, letter_heights, "(c)", fontsize=25)
    fig.text(0.7285, letter_heights, "(d)", fontsize=25)
    letter_heights = 0.4975
    fig.text(0.07125, letter_heights, "(e)", fontsize=25)
    fig.text(0.289, letter_heights, "(f)", fontsize=25)
    fig.text(0.509, letter_heights, "(g)", fontsize=25)
    fig.text(0.7285, letter_heights, "(h)", fontsize=25)
    os.system(f"mkdir -p {output_directory}")
    fig.savefig(f"{output_directory}/coulomb_poisson_global_twist_decomposition_with_emerg_field_"
                f"{output_file_string}.pdf", bbox_inches="tight")
    fig.clear()

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0))
    fig.tight_layout(w_pad=-12)
    plot_total_and_coulomb_spin_fields_of_single_row_fig(axes, lattice_coordinates, spin_phases, zero_temp_spin_phases)
    axes[2].quiver(lattice_coordinates[0], lattice_coordinates[1], np.cos(spin_wave_spin_phases),
                   np.sin(spin_wave_spin_phases), scale=25, width=0.005, color="r", pivot="mid")

    setup_axes_of_single_row_fig(axes)
    letter_heights = 0.975
    fig.text(0.09, letter_heights, "(a)", fontsize=25)
    fig.text(0.37, letter_heights, "(b)", fontsize=25)
    fig.text(0.6525, letter_heights, "(c)", fontsize=25)
    os.system(f"mkdir -p {output_directory}")
    fig.savefig(f"{output_directory}/coulomb_spin_wave_decomposition_{output_file_string}.pdf", bbox_inches="tight")
    fig.clear()

    fig, axes = plt.subplots(2, 3, figsize=(12.85, 8.0))
    fig.tight_layout(w_pad=-12)
    plot_total_and_coulomb_spin_fields_of_double_row_fig(axes, lattice_coordinates, spin_phases, zero_temp_spin_phases)
    axes[0, 2].quiver(lattice_coordinates[0], lattice_coordinates[1], np.cos(spin_wave_spin_phases),
                      np.sin(spin_wave_spin_phases), scale=25, width=0.005, color="r", pivot="mid")
    plot_total_and_coulomb_emergent_fields(axes, lattice_coordinates, electric_field_total, electric_field_irrotational,
                                           zeros)
    axes[1, 2].quiver(lattice_coordinates[0], lattice_coordinates[1] - 0.5, electric_field_auxiliary[0], zeros,
                      scale=20, width=0.002, color="k", pivot="mid")
    axes[1, 2].quiver(lattice_coordinates[0] - 0.5, lattice_coordinates[1], zeros, electric_field_auxiliary[1],
                      scale=20, width=0.002, color="k", pivot="mid")

    setup_axes_of_double_row_fig(axes)
    letter_heights = 0.99
    fig.text(0.09, letter_heights, "(a)", fontsize=25)
    fig.text(0.37, letter_heights, "(b)", fontsize=25)
    fig.text(0.6525, letter_heights, "(c)", fontsize=25)
    letter_heights = 0.4975
    fig.text(0.09, letter_heights, "(d)", fontsize=25)
    fig.text(0.37, letter_heights, "(e)", fontsize=25)
    fig.text(0.6525, letter_heights, "(f)", fontsize=25)
    os.system(f"mkdir -p {output_directory}")
    fig.savefig(f"{output_directory}/coulomb_spin_wave_decomposition_with_emerg_field_{output_file_string}.pdf",
                bbox_inches="tight")
    fig.clear()


def get_emergent_field_from_spin_phases(spin_phases, lattice_length):
    emergent_field = [[], []]
    for site in range(lattice_length ** 2):
        emergent_field[0].append(get_spin_difference(
            spin_phases[site], spin_phases[get_south_neighbour(site, lattice_length)]))
        emergent_field[1].append(-get_spin_difference(
            spin_phases[site], spin_phases[get_west_neighbour(site, lattice_length)]))
    return emergent_field


def get_spin_difference(spin_value_one, spin_value_two):
    return (spin_value_one - spin_value_two + math.pi) % (2.0 * math.pi) - math.pi


def plot_total_and_coulomb_spin_fields_of_single_row_fig(axes, lattice_coordinates, spin_phases, zero_temp_spin_phases):
    axes[0].quiver(lattice_coordinates[0], lattice_coordinates[1], np.cos(spin_phases), np.sin(spin_phases), scale=25,
                   width=0.005, color="r", pivot="mid")
    axes[1].quiver(lattice_coordinates[0], lattice_coordinates[1], np.cos(zero_temp_spin_phases),
                   np.sin(zero_temp_spin_phases), scale=25, width=0.005, color="r", pivot="mid")


def plot_total_and_coulomb_spin_fields_of_double_row_fig(axes, lattice_coordinates, spin_phases, zero_temp_spin_phases):
    axes[0, 0].quiver(lattice_coordinates[0], lattice_coordinates[1], np.cos(spin_phases), np.sin(spin_phases),
                      scale=25, width=0.005, color="r", pivot="mid")
    axes[0, 1].quiver(lattice_coordinates[0], lattice_coordinates[1], np.cos(zero_temp_spin_phases),
                      np.sin(zero_temp_spin_phases), scale=25, width=0.005, color="r", pivot="mid")


def plot_total_and_coulomb_emergent_fields(axes, lattice_coordinates, electric_field_total, electric_field_irrotational,
                                           zeros):
    axes[1, 0].quiver(lattice_coordinates[0], lattice_coordinates[1] - 0.5, electric_field_total[0], zeros,
                      scale=40, width=0.002, color="k", pivot="mid")
    axes[1, 0].quiver(lattice_coordinates[0] - 0.5, lattice_coordinates[1], zeros, electric_field_total[1],
                      scale=40, width=0.002, color="k", pivot="mid")
    axes[1, 1].quiver(lattice_coordinates[0], lattice_coordinates[1] - 0.5, electric_field_irrotational[0], zeros,
                      scale=40, width=0.002, color="k", pivot="mid")
    axes[1, 1].quiver(lattice_coordinates[0] - 0.5, lattice_coordinates[1], zeros, electric_field_irrotational[1],
                      scale=40, width=0.002, color="k", pivot="mid")


def setup_axes_of_single_row_fig(axes):
    [axis.axis("square") for axis in axes]
    [axis.xaxis.set_tick_params(labelbottom=False) for axis in axes]
    [axis.yaxis.set_tick_params(labelleft=False) for axis in axes]
    [axis.set_xticks([]) for axis in axes], [axis.set_yticks([]) for axis in axes]
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes]


def setup_axes_of_double_row_fig(axes):
    [axis.axis("square") for axis in axes.flat]
    [axis.xaxis.set_tick_params(labelbottom=False) for axis in axes.flat]
    [axis.yaxis.set_tick_params(labelleft=False) for axis in axes.flat]
    [axis.set_xticks([]) for axis in axes.flat], [axis.set_yticks([]) for axis in axes.flat]
    [axis.spines[spine].set_linewidth(3) for spine in ["top", "bottom", "left", "right"] for axis in axes.flat]


def get_quenched_spin_field(spin_field, integer_lattice_length=20, quench_with_xy_potential=True):
    for _ in range(2000):
        for i in range(integer_lattice_length ** 2):
            spin_perturbation = 0.01 * (np.random.uniform(0.0, 1.0) - 0.5)

            initial_spin_difference_1 = get_spin_difference(spin_field[i], spin_field[get_south_neighbour(i)])
            initial_spin_difference_2 = get_spin_difference(spin_field[get_west_neighbour(i)], spin_field[i])
            initial_spin_difference_3 = get_spin_difference(spin_field[get_north_neighbour(i)], spin_field[i])
            initial_spin_difference_4 = get_spin_difference(spin_field[i], spin_field[get_east_neighbour(i)])

            perturbed_spin_difference_1 = initial_spin_difference_1 + spin_perturbation
            perturbed_spin_difference_2 = initial_spin_difference_2 - spin_perturbation
            perturbed_spin_difference_3 = initial_spin_difference_3 - spin_perturbation
            perturbed_spin_difference_4 = initial_spin_difference_4 + spin_perturbation

            if ((abs(perturbed_spin_difference_1) < math.pi) and (abs(perturbed_spin_difference_2) < math.pi) and
                    (abs(perturbed_spin_difference_3) < math.pi) and (abs(perturbed_spin_difference_4) < math.pi)):
                # we may continue as candidate move does not alter the topological-defect configuration
                candidate_spin_value = (spin_field[i] + spin_perturbation) % (2.0 * math.pi)
                potential_difference = get_potential_difference(candidate_spin_value, spin_field, i,
                                                                quench_with_xy_potential)
                if potential_difference < 0.0:
                    spin_field[i] = candidate_spin_value
    return spin_field


def get_potential_difference(candidate_spin_value, spin_field, site_index, quench_with_xy_potential=True):
    if quench_with_xy_potential:
        return (- math.cos(spin_field[get_east_neighbour(site_index)] - candidate_spin_value) -
                math.cos(spin_field[get_north_neighbour(site_index)] - candidate_spin_value) -
                math.cos(candidate_spin_value - spin_field[get_west_neighbour(site_index)]) -
                math.cos(candidate_spin_value - spin_field[get_south_neighbour(site_index)]) +
                math.cos(spin_field[get_east_neighbour(site_index)] - spin_field[site_index]) +
                math.cos(spin_field[get_north_neighbour(site_index)] - spin_field[site_index]) +
                math.cos(spin_field[site_index] - spin_field[get_west_neighbour(site_index)]) +
                math.cos(spin_field[site_index] - spin_field[get_south_neighbour(site_index)]))
    # else, quench with HXY potential
    return 0.5 * (get_spin_difference(spin_field[get_east_neighbour(site_index)], candidate_spin_value) ** 2 +
                  get_spin_difference(spin_field[get_north_neighbour(site_index)], candidate_spin_value) ** 2 +
                  get_spin_difference(candidate_spin_value, spin_field[get_west_neighbour(site_index)]) ** 2 +
                  get_spin_difference(candidate_spin_value, spin_field[get_south_neighbour(site_index)]) ** 2 -
                  get_spin_difference(spin_field[get_east_neighbour(site_index)], spin_field[site_index]) ** 2 -
                  get_spin_difference(spin_field[get_north_neighbour(site_index)], spin_field[site_index]) ** 2 -
                  get_spin_difference(spin_field[site_index], spin_field[get_west_neighbour(site_index)]) ** 2 -
                  get_spin_difference(spin_field[site_index], spin_field[get_south_neighbour(site_index)]) ** 2)


def get_north_neighbour(i, integer_lattice_length=20):
    return i + ((int(i / integer_lattice_length) + 1) % integer_lattice_length -
                (int(i / integer_lattice_length)) % integer_lattice_length) * integer_lattice_length


def get_south_neighbour(i, integer_lattice_length=20):
    return i + ((int(i / integer_lattice_length) + integer_lattice_length - 1) % integer_lattice_length -
                (int(i / integer_lattice_length) + integer_lattice_length) % integer_lattice_length
                ) * integer_lattice_length


def get_east_neighbour(i, integer_lattice_length=20):
    return i + (i + 1) % integer_lattice_length - i % integer_lattice_length


def get_west_neighbour(i, integer_lattice_length=20):
    return i + (i - 1 + integer_lattice_length) % integer_lattice_length - (i + integer_lattice_length
                                                                            ) % integer_lattice_length


if __name__ == "__main__":
    main()
