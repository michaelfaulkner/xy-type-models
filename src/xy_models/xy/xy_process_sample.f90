subroutine process_sample(sample_index)
use variables
implicit none
integer :: i, sample_index, twist_relaxations(2), get_topological_sector_component, topological_sector(2)
double precision :: magnetic_norm_squared, potential, potential_cartesian_components(2), non_normalised_magnetisation(2)
double precision :: sum_of_1st_derivative_of_potential(2), sum_of_2nd_derivative_of_potential(2)
double precision :: sum_of_emergent_field(2), emergent_field_zero_mode(2), macro_josephson_current(2), inverse_vacuum_perm

if (measure_magnetisation) then
    non_normalised_magnetisation = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        non_normalised_magnetisation(1) = non_normalised_magnetisation(1) + cos(spin_field(i))
        non_normalised_magnetisation(2) = non_normalised_magnetisation(2) + sin(spin_field(i))
    end do
    if (print_samples) then
        write(30, 100) non_normalised_magnetisation(1), non_normalised_magnetisation(2)
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        ! nb, we divide by no_of_sites at this point to avoid infinities at large system size and long timescale
        magnetic_norm_squared = non_normalised_magnetisation(1) ** 2 / dfloat(no_of_sites ** 2) + &
                                non_normalised_magnetisation(2) ** 2 / dfloat(no_of_sites ** 2)
        magnetic_norm_sum = magnetic_norm_sum + magnetic_norm_squared ** 0.5
        magnetic_norm_squared_sum = magnetic_norm_squared_sum + magnetic_norm_squared
        magnetic_norm_quartic_sum = magnetic_norm_quartic_sum + magnetic_norm_squared ** 2
    end if
end if

if ((measure_helicity).or.(measure_potential).or.(measure_hot_twist_relaxations)) then
    potential_cartesian_components = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        potential_cartesian_components(1) = potential_cartesian_components(1) - &
                                                cos(spin_field(get_east_neighbour(i)) - spin_field(i))
        potential_cartesian_components(2) = potential_cartesian_components(2) - &
                                                cos(spin_field(get_north_neighbour(i)) - spin_field(i))
    end do
end if

if (measure_helicity) then
    sum_of_1st_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        sum_of_1st_derivative_of_potential(1) = sum_of_1st_derivative_of_potential(1) + &
                                                    sin(spin_field(get_east_neighbour(i)) - spin_field(i))
        sum_of_1st_derivative_of_potential(2) = sum_of_1st_derivative_of_potential(2) + &
                                                    sin(spin_field(get_north_neighbour(i)) - spin_field(i))
    end do
    sum_of_2nd_derivative_of_potential = (/ -potential_cartesian_components(1), -potential_cartesian_components(2) /)
    if (print_samples) then
        write(40, 100) sum_of_1st_derivative_of_potential(1), sum_of_1st_derivative_of_potential(2)
        write(50, 100) sum_of_2nd_derivative_of_potential(1), sum_of_2nd_derivative_of_potential(2)
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        ! nb, we divide by no_of_sites at this point to avoid infinities at large system size and long timescale
        inverse_vacuum_perm = 0.5d0 * (sum_of_2nd_derivative_of_potential(1) + sum_of_2nd_derivative_of_potential(2)) / &
                                        dfloat(no_of_sites)
        inverse_vacuum_perm_sum = inverse_vacuum_perm_sum + inverse_vacuum_perm
        inverse_vacuum_perm_squared_sum = inverse_vacuum_perm_squared_sum + inverse_vacuum_perm ** 2

        macro_josephson_current(1) = sum_of_1st_derivative_of_potential(1) / dfloat(no_of_sites)
        macro_josephson_current(2) = sum_of_1st_derivative_of_potential(2) / dfloat(no_of_sites)
        macro_josephson_current_sum(1) = macro_josephson_current_sum(1) + macro_josephson_current(1)
        macro_josephson_current_sum(2) = macro_josephson_current_sum(2) + macro_josephson_current(2)
        macro_josephson_current_squared_sum = macro_josephson_current_squared_sum + &
                                                 macro_josephson_current(1) ** 2 + macro_josephson_current(2) ** 2
        macro_josephson_current_quartic_sum = macro_josephson_current_quartic_sum + &
                                                (macro_josephson_current(1) ** 2 + macro_josephson_current(2) ** 2) ** 2
    end if
end if

if (measure_potential) then
    potential = potential_cartesian_components(1) + potential_cartesian_components(2)
    if (print_samples) then
        write(60, 200) potential
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        potential_sum = potential_sum + potential
        potential_squared_sum = potential_squared_sum + potential ** 2
        potential_quartic_sum = potential_quartic_sum + potential ** 4
    end if
end if

if (measure_hot_twist_relaxations) then
    twist_relaxations = (/ 0, 0 /)
    call twist_relaxations_calculation_x(twist_relaxations, 1, potential_cartesian_components) ! x direction; +ve twists
    call twist_relaxations_calculation_x(twist_relaxations, -1, potential_cartesian_components) ! x direction; -ve twists
    call twist_relaxations_calculation_y(twist_relaxations, 1, potential_cartesian_components) ! y direction; +ve twists
    call twist_relaxations_calculation_y(twist_relaxations, -1, potential_cartesian_components) ! y direction; -ve twists
    if (print_samples) then
        write(70, 300) twist_relaxations(1), twist_relaxations(2)
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        hot_twist_relaxations_sum(1) = hot_twist_relaxations_sum(1) + twist_relaxations(1)
        hot_twist_relaxations_sum(2) = hot_twist_relaxations_sum(2) + twist_relaxations(2)
        hot_twist_relaxations_squared_sum = hot_twist_relaxations_squared_sum + &
                                                twist_relaxations(1) ** 2 + twist_relaxations(2) ** 2
        hot_twist_relaxations_quartic_sum = hot_twist_relaxations_quartic_sum + &
                                                (twist_relaxations(1) ** 2 + twist_relaxations(2) ** 2) ** 2
    end if
end if

if ((measure_external_global_moves).and.(print_samples)) then
    write(80, 300) external_global_moves(1), external_global_moves(2)
end if

if (measure_twist_relaxations) then
    ! store spin field then slowly remove spin waves
    call store_spins
    do i = 1, 2000
        call xy_anneal_sweep
    end do

    ! now measure global twist-relaxation field (defined wrt annealed configuration)
    potential_cartesian_components = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        potential_cartesian_components(1) = potential_cartesian_components(1) - &
                                                cos(spin_field(get_east_neighbour(i)) - spin_field(i))
        potential_cartesian_components(2) = potential_cartesian_components(2) - &
                                                cos(spin_field(get_north_neighbour(i)) - spin_field(i))
    end do
    twist_relaxations = (/ 0, 0 /)
    call twist_relaxations_calculation_x(twist_relaxations, 1, potential_cartesian_components) ! x direction; +ve twists
    call twist_relaxations_calculation_x(twist_relaxations, -1, potential_cartesian_components) ! x direction; -ve twists
    call twist_relaxations_calculation_y(twist_relaxations, 1, potential_cartesian_components) ! y direction; +ve twists
    call twist_relaxations_calculation_y(twist_relaxations, -1, potential_cartesian_components) ! y direction; -ve twists
    if (print_samples) then
        write(90, 300) twist_relaxations(1), twist_relaxations(2)
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        twist_relaxations_sum(1) = twist_relaxations_sum(1) + twist_relaxations(1)
        twist_relaxations_sum(2) = twist_relaxations_sum(2) + twist_relaxations(2)
        twist_relaxations_squared_sum = twist_relaxations_squared_sum + &
                                            twist_relaxations(1) ** 2 + twist_relaxations(2) ** 2
        twist_relaxations_quartic_sum = twist_relaxations_quartic_sum + &
                                            (twist_relaxations(1) ** 2 + twist_relaxations(2) ** 2) ** 2
    end if

    ! restore spin field
    call get_spins
end if

if (measure_emergent_field) then
    call calculate_emergent_field
    sum_of_emergent_field = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        sum_of_emergent_field(1) = sum_of_emergent_field(1) + emergent_field(i, 1)
        sum_of_emergent_field(2) = sum_of_emergent_field(2) + emergent_field(i, 2)
    end do
    if (print_samples) then
        write(95, 100) sum_of_emergent_field(1), sum_of_emergent_field(2)
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        emergent_field_zero_mode(1) = sum_of_emergent_field(1) / dfloat(no_of_sites)
        emergent_field_zero_mode(2) = sum_of_emergent_field(2) / dfloat(no_of_sites)
        emergent_field_zero_mode_sum(1) = emergent_field_zero_mode_sum(1) + emergent_field_zero_mode(1)
        emergent_field_zero_mode_sum(2) = emergent_field_zero_mode_sum(2) + emergent_field_zero_mode(2)
        emergent_field_zero_mode_squared_sum = emergent_field_zero_mode_squared_sum + &
                                                 emergent_field_zero_mode(1) ** 2 + emergent_field_zero_mode(2) ** 2
        emergent_field_zero_mode_quartic_sum = emergent_field_zero_mode_quartic_sum + &
                                                (emergent_field_zero_mode(1) ** 2 + emergent_field_zero_mode(2) ** 2) ** 2

        topological_sector(1) = get_topological_sector_component(sum_of_emergent_field(1))
        topological_sector(2) = get_topological_sector_component(sum_of_emergent_field(2))
        topological_sector_sum(1) = topological_sector_sum(1) + topological_sector(1)
        topological_sector_sum(2) = topological_sector_sum(2) + topological_sector(2)
        topological_sector_squared_sum = topological_sector_squared_sum + &
                                            topological_sector(1) ** 2 + topological_sector(2) ** 2
        topological_sector_quartic_sum = topological_sector_quartic_sum + &
                                            (topological_sector(1) ** 2 + topological_sector(2) ** 2) ** 2
    end if
end if

100 format(ES25.14, ",", ES25.14)
200 format(ES25.14)
300 format(I10, ",", I10)

return
end subroutine process_sample


subroutine twist_relaxations_calculation_x(twist_relaxations, sign_of_twist, potential_cartesian_components)
use variables
implicit none
integer :: i, proposed_no_of_twists, twist_relaxations(2), sign_of_twist
double precision :: potential_difference, potential_of_twisted_field, potential_cartesian_components(2)
double precision :: potential_difference_w_previous_twist, potential_of_previous_twisted_field

potential_of_previous_twisted_field = potential_cartesian_components(1)
do proposed_no_of_twists = 1, integer_lattice_length - 1
    potential_of_twisted_field = 0.0d0
    do i = 1, no_of_sites
        potential_of_twisted_field = potential_of_twisted_field - &
                                        cos(spin_field(i) + dfloat(sign_of_twist) * dfloat(proposed_no_of_twists) * &
                                            two_pi / dfloat(integer_lattice_length) - spin_field(get_west_neighbour(i)))
    end do
    ! first check potential difference wrt previously twisted configuration - if 0.0 and sign_of_twist < 0, then exit;
    ! this accounts for energy degeneracy that can occur for a charge-neutral vortex pair separated by a distance L/2
    ! nb, sign_of_twist > 0 for y twists - the different choices align with the emergent-field convention
    potential_difference_w_previous_twist = potential_of_twisted_field - potential_of_previous_twisted_field
    if ((abs(potential_difference_w_previous_twist) < 1.0d-8).and.(sign_of_twist < 0)) then
       exit
    end if

    ! now check potential difference wrt original configuration
    potential_difference = potential_of_twisted_field - potential_cartesian_components(1)
    if ((potential_difference < 0.0d0).and.(abs(potential_difference) > 1.0d-12)) then
        twist_relaxations(1) = sign_of_twist * proposed_no_of_twists
    else
        exit
    end if
    potential_of_previous_twisted_field = potential_of_twisted_field
end do

return
end subroutine twist_relaxations_calculation_x


subroutine twist_relaxations_calculation_y(twist_relaxations, sign_of_twist, potential_cartesian_components)
use variables
implicit none
integer :: i, proposed_no_of_twists, twist_relaxations(2), sign_of_twist
double precision :: potential_difference, potential_of_twisted_field, potential_cartesian_components(2)
double precision :: potential_difference_w_previous_twist, potential_of_previous_twisted_field

potential_of_previous_twisted_field = potential_cartesian_components(2)
do proposed_no_of_twists = 1, integer_lattice_length - 1
    potential_of_twisted_field = 0.0d0
    do i = 1, no_of_sites
        potential_of_twisted_field = potential_of_twisted_field - &
                                        cos(spin_field(i) + dfloat(sign_of_twist) * dfloat(proposed_no_of_twists) * &
                                        two_pi / dfloat(integer_lattice_length) - spin_field(get_south_neighbour(i)))
    end do

    ! first check potential difference wrt previously twisted configuration - if 0.0 and sign_of_twist > 0, then exit;
    ! this accounts for energy degeneracy that can occur for a charge-neutral vortex pair separated by a distance L/2
    ! nb, sign_of_twist < 0 for x twists - the different choices align with the emergent-field convention
    potential_difference_w_previous_twist = potential_of_twisted_field - potential_of_previous_twisted_field
    if ((abs(potential_difference_w_previous_twist) < 1.0d-8).and.(sign_of_twist > 0)) then
       exit
    end if

    ! now check potential difference wrt original configuration
    potential_difference = potential_of_twisted_field - potential_cartesian_components(2)
    if ((potential_difference < 0.0d0).and.(abs(potential_difference) > 1.0d-12)) then
        twist_relaxations(2) = sign_of_twist * proposed_no_of_twists
    else
        exit
    end if
    potential_of_previous_twisted_field = potential_of_twisted_field
end do

return
end subroutine twist_relaxations_calculation_y


subroutine store_spins
use variables
implicit none
character(100) :: filename
integer :: spin_index, temperature_index

write(filename, '(A, "/stored_spins.csv")') trim(output_directory)
open(1000, file=filename)
do spin_index = 1, no_of_sites
    write(1000, 100) spin_field(spin_index)
end do
close(1000)

100 format(ES25.14)

return
end subroutine store_spins


subroutine get_spins
use variables
implicit none
character(100) :: filename
integer :: spin_index, temperature_index

write(filename, '(A, "/stored_spins.csv")') trim(output_directory)
open(1000, file=filename)
do spin_index = 1, no_of_sites
    read(1000, 100) spin_field(spin_index)
end do
close(1000)

100 format(ES25.14)

return
end subroutine get_spins
