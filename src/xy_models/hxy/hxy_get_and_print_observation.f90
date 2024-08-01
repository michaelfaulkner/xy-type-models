subroutine get_and_print_observation(observation_index)
use variables
implicit none
integer :: i, n, no_of_external_twists_to_minimise_potential(2), topological_sector(2), observation_index
double precision :: potential, sum_of_squared_emergent_field(2), non_normalised_magnetisation(2)
double precision :: sum_of_1st_derivative_of_potential(2), sum_of_2nd_derivative_of_potential(2)
double precision :: raw_magnetic_norm_squared, raw_inverse_vacuum_perm

if (measure_magnetisation) then
    non_normalised_magnetisation = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        non_normalised_magnetisation(1) = non_normalised_magnetisation(1) + cos(spin_field(i))
        non_normalised_magnetisation(2) = non_normalised_magnetisation(2) + sin(spin_field(i))
    end do
    if (print_samples) then
        write(30, 100) non_normalised_magnetisation(1), non_normalised_magnetisation(2)
    end if
    if (observation_index >= no_of_equilibration_sweeps) then
        ! magnetic_norm_squared = raw_magnetic_norm_squared / no_of_sites ** 2
        raw_magnetic_norm_squared = non_normalised_magnetisation(1) ** 2 + non_normalised_magnetisation(2) ** 2
        raw_magnetic_norm_sum = raw_magnetic_norm_sum + raw_magnetic_norm_squared ** 0.5
        raw_magnetic_norm_squared_sum = raw_magnetic_norm_squared_sum + raw_magnetic_norm_squared
        raw_magnetic_norm_quartic_sum = raw_magnetic_norm_quartic_sum + raw_magnetic_norm_squared ** 2
    end if
end if

if (measure_helicity) then
    sum_of_1st_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
    sum_of_2nd_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        sum_of_1st_derivative_of_potential(1) = sum_of_1st_derivative_of_potential(1) + emergent_field(i, 2)
        sum_of_1st_derivative_of_potential(2) = sum_of_1st_derivative_of_potential(2) + emergent_field(i, 1)
        do n = 1, vacuum_permittivity_sum_cutoff
            sum_of_2nd_derivative_of_potential(1) = sum_of_2nd_derivative_of_potential(1) + &
                                                            (-1.0d0) ** (n + 1) * cos(dfloat(n) * emergent_field(i, 2))
            sum_of_2nd_derivative_of_potential(2) = sum_of_2nd_derivative_of_potential(2) + &
                                                            (-1.0d0) ** (n + 1) * cos(dfloat(n) * emergent_field(i, 1))
        end do
    end do
    sum_of_2nd_derivative_of_potential(1) = 2.0d0 * sum_of_2nd_derivative_of_potential(1)
    sum_of_2nd_derivative_of_potential(2) = 2.0d0 * sum_of_2nd_derivative_of_potential(2)
    if (print_samples) then
        write(40, 100) sum_of_1st_derivative_of_potential(1), sum_of_1st_derivative_of_potential(2)
        write(50, 100) sum_of_2nd_derivative_of_potential(1), sum_of_2nd_derivative_of_potential(2)
    end if
    if (observation_index >= no_of_equilibration_sweeps) then
        ! inverse_vacuum_perm = raw_inverse_vacuum_perm / 2
        raw_inverse_vacuum_perm = sum_of_2nd_derivative_of_potential(1) + sum_of_2nd_derivative_of_potential(2)
        raw_inverse_vacuum_perm_sum = raw_inverse_vacuum_perm_sum + raw_inverse_vacuum_perm
        raw_inverse_vacuum_perm_squared_sum = raw_inverse_vacuum_perm_squared_sum + raw_inverse_vacuum_perm ** 2
        ! macro_josephson_current = raw_macro_josephson_current / no_of_sites
        raw_macro_josephson_current_sum(1) = raw_macro_josephson_current_sum(1) + sum_of_1st_derivative_of_potential(1)
        raw_macro_josephson_current_sum(2) = raw_macro_josephson_current_sum(2) + sum_of_1st_derivative_of_potential(2)
        raw_macro_josephson_current_squared_sum = raw_macro_josephson_current_squared_sum + &
                                sum_of_1st_derivative_of_potential(1) ** 2 + sum_of_1st_derivative_of_potential(2) ** 2
        raw_macro_josephson_current_quartic_sum = raw_macro_josephson_current_quartic_sum + &
                        (sum_of_1st_derivative_of_potential(1) ** 2 + sum_of_1st_derivative_of_potential(2) ** 2) ** 2
        ! w_{x / y} = floor((sum_r E_{x / y}(r) + pi L) / (2pi L)), where w \in Z^2 is the topological_sector
        topological_sector(1) = floor(sum_of_1st_derivative_of_potential(2) / two_pi / integer_lattice_length + 0.5d0)
        topological_sector(2) = floor(-sum_of_1st_derivative_of_potential(1) / two_pi / integer_lattice_length + 0.5d0)
        topological_sector_sum(1) = topological_sector_sum(1) + topological_sector(1)
        topological_sector_sum(2) = topological_sector_sum(2) + topological_sector(2)
        topological_sector_squared_sum = topological_sector_squared_sum + &
                                            topological_sector(1) ** 2 + topological_sector(2) ** 2
        topological_sector_quartic_sum = topological_sector_quartic_sum + &
                                            (topological_sector(1) ** 2 + topological_sector(2) ** 2) ** 2
    end if
end if

if ((measure_potential).or.(measure_potential_minimising_twists)) then
    sum_of_squared_emergent_field = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        sum_of_squared_emergent_field(1) = sum_of_squared_emergent_field(1) + emergent_field(i, 1) ** 2
        sum_of_squared_emergent_field(2) = sum_of_squared_emergent_field(2) + emergent_field(i, 2) ** 2
    end do
end if

if (measure_potential) then
    potential = 0.5d0 * (sum_of_squared_emergent_field(1) + sum_of_squared_emergent_field(2))
    if (print_samples) then
        write(60, 200) potential
    end if
    if (observation_index >= no_of_equilibration_sweeps) then
        potential_sum = potential_sum + potential
        potential_squared_sum = potential_squared_sum + potential ** 2
    end if
end if

if (measure_potential_minimising_twists) then
    no_of_external_twists_to_minimise_potential = (/ 0, 0 /)
    call potential_minimising_twists_calculation_x(no_of_external_twists_to_minimise_potential, 1, &
                                                        sum_of_squared_emergent_field) ! x direction; positive twists
    call potential_minimising_twists_calculation_x(no_of_external_twists_to_minimise_potential, -1, &
                                                        sum_of_squared_emergent_field) ! x direction; negative twists
    call potential_minimising_twists_calculation_y(no_of_external_twists_to_minimise_potential, 1, &
                                                        sum_of_squared_emergent_field) ! y direction; positive twists
    call potential_minimising_twists_calculation_y(no_of_external_twists_to_minimise_potential, -1, &
                                                        sum_of_squared_emergent_field) ! y direction; negative twists
    if (print_samples) then
        write(70, 300) no_of_external_twists_to_minimise_potential(1), no_of_external_twists_to_minimise_potential(2)
    end if
    if (observation_index >= no_of_equilibration_sweeps) then
        potential_minimising_twists_sum(1) = potential_minimising_twists_sum(1) + &
                                                no_of_external_twists_to_minimise_potential(1)
        potential_minimising_twists_sum(2) = potential_minimising_twists_sum(2) + &
                                                no_of_external_twists_to_minimise_potential(2)
        potential_minimising_twists_squared_sum = potential_minimising_twists_squared_sum + &
                no_of_external_twists_to_minimise_potential(1) ** 2 + no_of_external_twists_to_minimise_potential(2) ** 2
        potential_minimising_twists_quartic_sum = potential_minimising_twists_quartic_sum + &
                (no_of_external_twists_to_minimise_potential(1) ** 2 + no_of_external_twists_to_minimise_potential(2) ** 2) ** 2
    end if
end if

if ((measure_external_global_moves).and.(print_samples)) then
    write(80, 300) external_global_moves(1), external_global_moves(2)
end if

100 format(ES25.14, ",", ES25.14)
200 format(ES25.14)
300 format(I10, ",", I10)

return
end subroutine get_and_print_observation


subroutine potential_minimising_twists_calculation_x(no_of_external_twists_to_minimise_potential, sign_of_twist, &
                                                        sum_of_squared_emergent_field)
use variables
implicit none
integer :: i, proposed_no_of_twists, no_of_external_twists_to_minimise_potential(2), sign_of_twist
double precision :: potential_difference, get_spin_difference, sum_of_squared_twisted_emergent_field
double precision :: sum_of_squared_emergent_field(2)

proposed_no_of_twists = 1
do
    sum_of_squared_twisted_emergent_field = 0.0d0
    do i = 1, no_of_sites
        sum_of_squared_twisted_emergent_field = sum_of_squared_twisted_emergent_field + &
                                                    get_spin_difference(spin_field(i) + dfloat(sign_of_twist) * &
                                                                        dfloat(proposed_no_of_twists) * two_pi / &
                                                                        dfloat(integer_lattice_length), &
                                                                        spin_field(get_west_neighbour(i))) ** 2
    end do
    potential_difference = 0.5d0 * (sum_of_squared_twisted_emergent_field - sum_of_squared_emergent_field(2))
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential(1) = sign_of_twist * proposed_no_of_twists
        proposed_no_of_twists = proposed_no_of_twists + 1
    else
        exit
    end if
end do

return
end subroutine potential_minimising_twists_calculation_x


subroutine potential_minimising_twists_calculation_y(no_of_external_twists_to_minimise_potential, sign_of_twist, &
                                                        sum_of_squared_emergent_field)
use variables
implicit none
integer :: i, proposed_no_of_twists, no_of_external_twists_to_minimise_potential(2), sign_of_twist
double precision :: potential_difference, get_spin_difference, sum_of_squared_twisted_emergent_field
double precision :: sum_of_squared_emergent_field(2)

proposed_no_of_twists = 1
do
    sum_of_squared_twisted_emergent_field = 0.0d0
    do i = 1, no_of_sites
        sum_of_squared_twisted_emergent_field = sum_of_squared_twisted_emergent_field + &
                                                    get_spin_difference(spin_field(i) + dfloat(sign_of_twist) * &
                                                                        dfloat(proposed_no_of_twists) * two_pi / &
                                                                        dfloat(integer_lattice_length), &
                                                                        spin_field(get_south_neighbour(i))) ** 2
    end do
    potential_difference = 0.5d0 * (sum_of_squared_twisted_emergent_field - sum_of_squared_emergent_field(1))
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential(2) = sign_of_twist * proposed_no_of_twists
        proposed_no_of_twists = proposed_no_of_twists + 1
    else
        exit
    end if
end do

return
end subroutine potential_minimising_twists_calculation_y
