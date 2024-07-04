subroutine get_and_print_observation
use variables
implicit none
integer :: i, no_of_global_twists_to_minimise_potential(2)
double precision :: potential, potential_cartesian_components(2), non_normalised_magnetisation(2)
double precision :: sum_of_1st_derivative_of_potential(2), sum_of_2nd_derivative_of_potential(2)
double precision :: sum_of_emergent_field(2)

if (measure_magnetisation) then
    non_normalised_magnetisation = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        non_normalised_magnetisation(1) = non_normalised_magnetisation(1) + cos(spin_field(i))
        non_normalised_magnetisation(2) = non_normalised_magnetisation(2) + sin(spin_field(i))
    end do
    write(30, 100) non_normalised_magnetisation(1), non_normalised_magnetisation(2)
end if

if ((measure_helicity).or.(measure_potential).or.(measure_potential_minimising_twists)) then
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
    write(40, 100) sum_of_1st_derivative_of_potential(1), sum_of_1st_derivative_of_potential(2)
    write(50, 100) sum_of_2nd_derivative_of_potential(1), sum_of_2nd_derivative_of_potential(2)
end if

if (measure_potential) then
    potential = potential_cartesian_components(1) + potential_cartesian_components(2)
    write(60, 200) potential
end if

if (measure_potential_minimising_twists) then
    no_of_global_twists_to_minimise_potential = (/ 0, 0 /)
    call potential_minimising_twists_calculation_x(no_of_global_twists_to_minimise_potential, 1, &
                                                    potential_cartesian_components) ! x direction; positive twists
    call potential_minimising_twists_calculation_x(no_of_global_twists_to_minimise_potential, -1, &
                                                    potential_cartesian_components) ! x direction; negative twists
    call potential_minimising_twists_calculation_y(no_of_global_twists_to_minimise_potential, 1, &
                                                    potential_cartesian_components) ! y direction; positive twists
    call potential_minimising_twists_calculation_y(no_of_global_twists_to_minimise_potential, -1, &
                                                    potential_cartesian_components) ! y direction; negative twists
    write(70, 300) no_of_global_twists_to_minimise_potential(1), no_of_global_twists_to_minimise_potential(2)
end if

if (measure_external_global_moves) then
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
    no_of_global_twists_to_minimise_potential = (/ 0, 0 /)
    call potential_minimising_twists_calculation_x(no_of_global_twists_to_minimise_potential, 1, &
                                                    potential_cartesian_components) ! x direction; positive twists
    call potential_minimising_twists_calculation_x(no_of_global_twists_to_minimise_potential, -1, &
                                                    potential_cartesian_components) ! x direction; negative twists
    call potential_minimising_twists_calculation_y(no_of_global_twists_to_minimise_potential, 1, &
                                                    potential_cartesian_components) ! y direction; positive twists
    call potential_minimising_twists_calculation_y(no_of_global_twists_to_minimise_potential, -1, &
                                                    potential_cartesian_components) ! y direction; negative twists
    write(90, 300) no_of_global_twists_to_minimise_potential(1), no_of_global_twists_to_minimise_potential(2)

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
    write(95, 100) sum_of_emergent_field(1), sum_of_emergent_field(2)
end if

100 format(ES25.14, ",", ES25.14)
200 format(ES25.14)
300 format(I10, ",", I10)

return
end subroutine get_and_print_observation


subroutine potential_minimising_twists_calculation_x(no_of_global_twists_to_minimise_potential, sign_of_twist, &
                                                        potential_cartesian_components)
use variables
implicit none
integer :: i, proposed_no_of_twists, no_of_global_twists_to_minimise_potential(2), sign_of_twist
double precision :: potential_difference, potential_of_twisted_field, potential_cartesian_components(2)

proposed_no_of_twists = 1
do
    potential_of_twisted_field = 0.0d0
    do i = 1, no_of_sites
        potential_of_twisted_field = potential_of_twisted_field - &
                                        cos(spin_field(i) + dfloat(sign_of_twist) * dfloat(proposed_no_of_twists) * &
                                            two_pi / dfloat(integer_lattice_length) - spin_field(get_west_neighbour(i)))
    end do
    potential_difference = potential_of_twisted_field - potential_cartesian_components(1)
    ! nb, we need abs(potential_difference) > 1.0d-12 below for cases où proposed_no_of_twists = integer_lattice_length
    ! in such cases, potential_difference is exactly zero leading to floating-point errors - this can be a problem at
    ! L = 4 as proposed_no_of_twists = 4 is reachable given one of {cos(n 2pi / 4), sin(n 2pi / 4)} is zero for all n
    if ((potential_difference < 0.0d0).and.(abs(potential_difference) > 1.0d-12)) then
        no_of_global_twists_to_minimise_potential(1) = sign_of_twist * proposed_no_of_twists
        proposed_no_of_twists = proposed_no_of_twists + 1
    else
        exit
    end if
end do

return
end subroutine potential_minimising_twists_calculation_x


subroutine potential_minimising_twists_calculation_y(no_of_global_twists_to_minimise_potential, sign_of_twist, &
                                                        potential_cartesian_components)
use variables
implicit none
integer :: i, proposed_no_of_twists, no_of_global_twists_to_minimise_potential(2), sign_of_twist
double precision :: potential_difference, potential_of_twisted_field, potential_cartesian_components(2)

proposed_no_of_twists = 1
do
    potential_of_twisted_field = 0.0d0
    do i = 1, no_of_sites
        potential_of_twisted_field = potential_of_twisted_field - &
                                        cos(spin_field(i) + dfloat(sign_of_twist) * dfloat(proposed_no_of_twists) * &
                                        two_pi / dfloat(integer_lattice_length) - spin_field(get_south_neighbour(i)))
    end do
    potential_difference = potential_of_twisted_field - potential_cartesian_components(2)
    ! nb, we need abs(potential_difference) > 1.0d-12 below for cases où proposed_no_of_twists = integer_lattice_length
    ! in such cases, potential_difference is exactly zero leading to floating-point errors - this can be a problem at
    ! L = 4 as proposed_no_of_twists = 4 is reachable given one of {cos(n 2pi / 4), sin(n 2pi / 4)} is zero for all n
    if ((potential_difference < 0.0d0).and.(abs(potential_difference) > 1.0d-12)) then
        no_of_global_twists_to_minimise_potential(2) = sign_of_twist * proposed_no_of_twists
        proposed_no_of_twists = proposed_no_of_twists + 1
    else
        exit
    end if
end do

return
end subroutine potential_minimising_twists_calculation_y


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
