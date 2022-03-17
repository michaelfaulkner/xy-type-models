subroutine get_and_print_observation
use variables
implicit none
integer :: i, no_of_external_twists_to_minimise_potential(2)
double precision :: potential, potential_cartesian_components(2), non_normalised_magnetisation(2)
double precision :: sum_of_1st_derivative_of_potential(2), sum_of_2nd_derivative_of_potential(2)

no_of_external_twists_to_minimise_potential = (/ 0, 0 /)
non_normalised_magnetisation = (/ 0.0d0, 0.0d0 /)
sum_of_1st_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
sum_of_2nd_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
potential_cartesian_components = (/ 0.0d0, 0.0d0 /)

do i = 1, no_of_sites
    non_normalised_magnetisation(1) = non_normalised_magnetisation(1) + cos(spin_field(i))
    non_normalised_magnetisation(2) = non_normalised_magnetisation(2) + sin(spin_field(i))
    sum_of_1st_derivative_of_potential(1) = sum_of_1st_derivative_of_potential(1) + &
                                                sin(spin_field(get_east_neighbour(i)) - spin_field(i))
    sum_of_1st_derivative_of_potential(2) = sum_of_1st_derivative_of_potential(2) + &
                                                sin(spin_field(get_north_neighbour(i)) - spin_field(i))
    sum_of_2nd_derivative_of_potential(1) = sum_of_2nd_derivative_of_potential(1) + &
                                                cos(spin_field(get_east_neighbour(i)) - spin_field(i))
    sum_of_2nd_derivative_of_potential(2) = sum_of_2nd_derivative_of_potential(2) + &
                                                cos(spin_field(get_north_neighbour(i)) - spin_field(i))
    potential_cartesian_components(1) = potential_cartesian_components(1) - &
                                            cos(spin_field(get_east_neighbour(i)) - spin_field(i))
    potential_cartesian_components(2) = potential_cartesian_components(2) - &
                                            cos(spin_field(get_north_neighbour(i)) - spin_field(i))
end do
potential = potential_cartesian_components(1) + potential_cartesian_components(2)

write(20, 100) potential
if (use_external_global_moves) then
    write(30, 200) external_global_move(1), external_global_move(2)
end if
write(40, 300) non_normalised_magnetisation(1), non_normalised_magnetisation(2)
write(50, 300) sum_of_1st_derivative_of_potential(1), sum_of_1st_derivative_of_potential(2)
write(60, 300) sum_of_2nd_derivative_of_potential(1), sum_of_2nd_derivative_of_potential(2)

if (calculate_potential_minimising_twists) then
    call potential_minimising_twists_calculation_x(no_of_external_twists_to_minimise_potential, 1, &
                                                    potential_cartesian_components) ! x direction; positive twists
    call potential_minimising_twists_calculation_x(no_of_external_twists_to_minimise_potential, -1, &
                                                    potential_cartesian_components) ! x direction; negative twists
    call potential_minimising_twists_calculation_y(no_of_external_twists_to_minimise_potential, 1, &
                                                    potential_cartesian_components) ! y direction; positive twists
    call potential_minimising_twists_calculation_y(no_of_external_twists_to_minimise_potential, -1, &
                                                    potential_cartesian_components) ! y direction; negative twists
    write(70, 200) no_of_external_twists_to_minimise_potential(1), no_of_external_twists_to_minimise_potential(2)
end if

100 format(ES25.14)
200 format(I10, ",", I10)
300 format(ES25.14, ",", ES25.14)

return
end subroutine get_and_print_observation


subroutine potential_minimising_twists_calculation_x(no_of_external_twists_to_minimise_potential, sign_of_twist, &
                                                        potential_cartesian_components)
use variables
implicit none
integer :: i, proposed_no_of_twists, no_of_external_twists_to_minimise_potential(2), sign_of_twist
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
                                                        potential_cartesian_components)
use variables
implicit none
integer :: i, proposed_no_of_twists, no_of_external_twists_to_minimise_potential(2), sign_of_twist
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
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential(2) = sign_of_twist * proposed_no_of_twists
        proposed_no_of_twists = proposed_no_of_twists + 1
    else
        exit
    end if
end do

return
end subroutine potential_minimising_twists_calculation_y
