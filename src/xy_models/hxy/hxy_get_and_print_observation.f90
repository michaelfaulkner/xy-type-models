subroutine get_and_print_observation
use variables
implicit none
integer :: i, n
double precision :: potential, non_normalised_magnetisation(2)
double precision :: sum_of_1st_derivative_of_potential(2), sum_of_2nd_derivative_of_potential(2)

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
                                                                    (-1.0d0) ** (n + 1) * cos(n * emergent_field(i, 2))
        sum_of_2nd_derivative_of_potential(2) = sum_of_2nd_derivative_of_potential(2) + &
                                                                    (-1.0d0) ** (n + 1) * cos(n * emergent_field(i, 1))
    end do
    sum_of_squared_emergent_field(1) = sum_of_squared_emergent_field(1) + emergent_field(i, 1) * emergent_field(i, 1)
    sum_of_squared_emergent_field(2) = sum_of_squared_emergent_field(2) + emergent_field(i, 2) * emergent_field(i, 2)
end do
sum_of_2nd_derivative_of_potential(1) = 2.0d0 * sum_of_2nd_derivative_of_potential(1)
sum_of_2nd_derivative_of_potential(2) = 2.0d0 * sum_of_2nd_derivative_of_potential(2)
potential = 0.5d0 * (sum_of_squared_emergent_field(1) + sum_of_squared_emergent_field(2))

if (calculate_external_minimising_twist_field) then
    call external_minimising_twist_field_calculation
end if

write(20, 100) potential, non_normalised_magnetisation(1), non_normalised_magnetisation(2), &
                sum_of_1st_derivative_of_potential(1), sum_of_1st_derivative_of_potential(2), &
                sum_of_2nd_derivative_of_potential(1), sum_of_2nd_derivative_of_potential(2), &
                external_global_move(1), external_global_move(2), &
                no_of_external_twists_to_minimise_potential(1), no_of_external_twists_to_minimise_potential(2)

100 format(ES30.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", I29, &
            ",", I29, ",", I29, ",", I29)

return
end subroutine get_and_print_observation


subroutine external_minimising_twist_field_calculation
use variables
implicit none
integer :: i, n
double precision :: potential_difference, get_spin_difference, current_sum_of_squared_emergent_field
double precision :: twisted_sum_of_squared_emergent_field

! x direction; positive twist
n = 1
do
    current_sum_of_squared_emergent_field = sum_of_squared_emergent_field(2)
    twisted_sum_of_squared_emergent_field = 0.0d0
    do i = 1, no_of_sites
        twisted_sum_of_squared_emergent_field = twisted_sum_of_squared_emergent_field &
                            + get_spin_difference(spin_field(i) + dfloat(n) * two_pi / dfloat(integer_lattice_length), &
                                                  spin_field(get_west_neighbour(i))) ** 2
    end do
    potential_difference = 0.5d0 * (twisted_sum_of_squared_emergent_field - current_sum_of_squared_emergent_field)
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential(1) = n
        n = n + 1
    else
        exit
    end if
end do

! x direction; negative twist
n = 1
do
    current_sum_of_squared_emergent_field = sum_of_squared_emergent_field(2)
    twisted_sum_of_squared_emergent_field = 0.0d0
    do i = 1, no_of_sites
        twisted_sum_of_squared_emergent_field = twisted_sum_of_squared_emergent_field &
                            + get_spin_difference(spin_field(i) - dfloat(n) * two_pi / dfloat(integer_lattice_length), &
                                                  spin_field(get_west_neighbour(i)))** 2
    end do
    potential_difference = 0.5d0 * (twisted_sum_of_squared_emergent_field - current_sum_of_squared_emergent_field)
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential(1) = n
        n = n + 1
    else
        exit
    end if
end do

! y direction; positive twist
n = 1
do
    current_sum_of_squared_emergent_field = sum_of_squared_emergent_field(1)
    twisted_sum_of_squared_emergent_field = 0.0d0
    do i = 1, no_of_sites
        twisted_sum_of_squared_emergent_field = twisted_sum_of_squared_emergent_field &
                            + get_spin_difference(spin_field(i) + dfloat(n) * two_pi / dfloat(integer_lattice_length), &
                                                  spin_field(get_south_neighbour(i))) ** 2
    end do
    potential_difference = 0.5d0 * (twisted_sum_of_squared_emergent_field - current_sum_of_squared_emergent_field)
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential(2) = n
        n = n + 1
    else
        exit
    end if
end do

! y direction; negative twist
n = 1
do
    current_sum_of_squared_emergent_field = sum_of_squared_emergent_field(1)
    twisted_sum_of_squared_emergent_field = 0.0d0
    do i = 1, no_of_sites
        twisted_sum_of_squared_emergent_field = twisted_sum_of_squared_emergent_field &
                            + get_spin_difference(spin_field(i) - dfloat(n) * two_pi / dfloat(integer_lattice_length), &
                                                  spin_field(get_south_neighbour(i))) ** 2
    end do
    potential_difference = 0.5d0 * (twisted_sum_of_squared_emergent_field - current_sum_of_squared_emergent_field)
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential(2) = n
        n = n + 1
    else
        exit
    end if
end do

return
end subroutine
