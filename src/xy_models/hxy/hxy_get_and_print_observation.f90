subroutine get_and_print_observation
use variables
implicit none
integer :: i, n, no_of_external_twists_to_minimise_potential(2)
double precision :: potential, sum_of_squared_emergent_field(2), non_normalised_magnetisation(2)
double precision :: sum_of_1st_derivative_of_potential(2), sum_of_2nd_derivative_of_potential(2)

no_of_external_twists_to_minimise_potential = (/ 0, 0 /)
non_normalised_magnetisation = (/ 0.0d0, 0.0d0 /)
sum_of_1st_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
sum_of_2nd_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
sum_of_squared_emergent_field = (/ 0.0d0, 0.0d0 /)

do i = 1, no_of_sites
    non_normalised_magnetisation(1) = non_normalised_magnetisation(1) + cos(spin_field(i))
    non_normalised_magnetisation(2) = non_normalised_magnetisation(2) + sin(spin_field(i))
    sum_of_1st_derivative_of_potential(1) = sum_of_1st_derivative_of_potential(1) + emergent_field(i, 2)
    sum_of_1st_derivative_of_potential(2) = sum_of_1st_derivative_of_potential(2) + emergent_field(i, 1)
    do n = 1, vacuum_permittivity_sum_cutoff
        sum_of_2nd_derivative_of_potential(1) = sum_of_2nd_derivative_of_potential(1) + &
                                                            (-1.0d0) ** (n + 1) * cos(dfloat(n) * emergent_field(i, 2))
        sum_of_2nd_derivative_of_potential(2) = sum_of_2nd_derivative_of_potential(2) + &
                                                            (-1.0d0) ** (n + 1) * cos(dfloat(n) * emergent_field(i, 1))
    end do
    sum_of_squared_emergent_field(1) = sum_of_squared_emergent_field(1) + emergent_field(i, 1) ** 2
    sum_of_squared_emergent_field(2) = sum_of_squared_emergent_field(2) + emergent_field(i, 2) ** 2
end do
sum_of_2nd_derivative_of_potential(1) = 2.0d0 * sum_of_2nd_derivative_of_potential(1)
sum_of_2nd_derivative_of_potential(2) = 2.0d0 * sum_of_2nd_derivative_of_potential(2)
potential = 0.5d0 * (sum_of_squared_emergent_field(1) + sum_of_squared_emergent_field(2))

write(20, 100) potential
write(30, 200) non_normalised_magnetisation(1), non_normalised_magnetisation(2)
write(40, 200) sum_of_1st_derivative_of_potential(1), sum_of_1st_derivative_of_potential(2)
write(50, 200) sum_of_2nd_derivative_of_potential(1), sum_of_2nd_derivative_of_potential(2)
if (use_external_global_moves) then
    write(60, 300) external_global_moves(1), external_global_moves(2)
end if

if (calculate_potential_minimising_twists) then
    call potential_minimising_twists_calculation_x(no_of_external_twists_to_minimise_potential, 1, &
                                                        sum_of_squared_emergent_field) ! x direction; positive twists
    call potential_minimising_twists_calculation_x(no_of_external_twists_to_minimise_potential, -1, &
                                                        sum_of_squared_emergent_field) ! x direction; negative twists
    call potential_minimising_twists_calculation_y(no_of_external_twists_to_minimise_potential, 1, &
                                                        sum_of_squared_emergent_field) ! y direction; positive twists
    call potential_minimising_twists_calculation_y(no_of_external_twists_to_minimise_potential, -1, &
                                                        sum_of_squared_emergent_field) ! y direction; negative twists
    write(70, 300) no_of_external_twists_to_minimise_potential(1), no_of_external_twists_to_minimise_potential(2)
end if

100 format(ES25.14)
200 format(ES25.14, ",", ES25.14)
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
