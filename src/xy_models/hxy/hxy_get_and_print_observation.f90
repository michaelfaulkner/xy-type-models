subroutine get_and_print_observation
use variables
implicit none
integer :: i, n
double precision :: potential
double precision :: magnetisation(2), sum_of_1st_derivative_of_potential(2), sum_of_2nd_derivative_of_potential(2)

magnetisation = (/ 0.0d0, 0.0d0 /)
sum_of_1st_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
sum_of_2nd_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
sum_of_squared_electric_field_x = 0.0d0
sum_of_squared_electric_field_y = 0.0d0

do i = 1, no_of_sites
    magnetisation(1) = magnetisation(1) + cos(spin_field(i))
    magnetisation(2) = magnetisation(2) + sin(spin_field(i))
    sum_of_1st_derivative_of_potential(1) = sum_of_1st_derivative_of_potential(1) + emergent_field_y(i)
    sum_of_1st_derivative_of_potential(2) = sum_of_1st_derivative_of_potential(2) + emergent_field_x(i)
    do n = 1, vacuum_permittivity_sum_cutoff
        sum_of_2nd_derivative_of_potential(1) = sum_of_2nd_derivative_of_potential(1) + &
                                                                    (-1.0d0) ** (n + 1) * cos(n * emergent_field_y(i))
        sum_of_2nd_derivative_of_potential(2) = sum_of_2nd_derivative_of_potential(2) + &
                                                                    (-1.0d0) ** (n + 1) * cos(n * emergent_field_x(i))
    end do
    sum_of_squared_electric_field_x = sum_of_squared_electric_field_x + emergent_field_x(i) * emergent_field_x(i)
    sum_of_squared_electric_field_y = sum_of_squared_electric_field_y + emergent_field_y(i) * emergent_field_y(i)
end do
sum_of_2nd_derivative_of_potential(1) = 2.0d0 * sum_of_2nd_derivative_of_potential(1)
sum_of_2nd_derivative_of_potential(2) = 2.0d0 * sum_of_2nd_derivative_of_potential(2)
potential = 0.5d0 * (sum_of_squared_electric_field_x + sum_of_squared_electric_field_y)

if (calculate_external_minimising_twist_field) then
    call external_minimising_twist_field_calculation
end if
  
write(20, 100) potential, magnetisation(1), magnetisation(2), &
                sum_of_1st_derivative_of_potential(1), sum_of_1st_derivative_of_potential(2), &
                sum_of_2nd_derivative_of_potential(1), sum_of_2nd_derivative_of_potential(2), &
                no_of_external_twists_to_minimise_potential_x, no_of_external_twists_to_minimise_potential_y

100 format(ES30.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", I29, &
            ",", I29)

return
end subroutine get_and_print_observation


subroutine external_minimising_twist_field_calculation
use variables
implicit none
integer i, n
double precision current_sum_of_squared_electric_field_x, current_sum_of_squared_electric_field_y
double precision twisted_sum_of_squared_electric_field_x, twisted_sum_of_squared_electric_field_y, potential_difference

! x direction; positive twist
n = 1
do
    current_sum_of_squared_electric_field_y = sum_of_squared_electric_field_y
    twisted_sum_of_squared_electric_field_y = 0.0d0
    do i = 1, no_of_sites
        twisted_sum_of_squared_electric_field_y = twisted_sum_of_squared_electric_field_y &
                + (modulo(spin_field(i) - spin_field(get_west_neighbour(i)) &
                        + n * twopi / integer_lattice_length + pi, twopi) - pi) ** 2
    end do
    potential_difference = 0.5d0 * (twisted_sum_of_squared_electric_field_y - current_sum_of_squared_electric_field_y)
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential_x = n
        n = n + 1
    else
        exit
    end if
end do

! x direction; negative twist
n = 1
do
    current_sum_of_squared_electric_field_y = sum_of_squared_electric_field_y
    twisted_sum_of_squared_electric_field_y = 0.0d0
    do i = 1, no_of_sites
        twisted_sum_of_squared_electric_field_y = twisted_sum_of_squared_electric_field_y &
                + (modulo(spin_field(i) - spin_field(get_west_neighbour(i)) &
                        - n * twopi / integer_lattice_length + pi, twopi) - pi) ** 2
    end do
    potential_difference = 0.5d0 * (twisted_sum_of_squared_electric_field_y - current_sum_of_squared_electric_field_y)
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential_x = n
        n = n + 1
    else
        exit
    end if
end do

! y direction; positive twist
n = 1
do
    current_sum_of_squared_electric_field_x = sum_of_squared_electric_field_x
    twisted_sum_of_squared_electric_field_x = 0.0d0
    do i = 1, no_of_sites
        twisted_sum_of_squared_electric_field_x = twisted_sum_of_squared_electric_field_x &
                + (modulo(spin_field(i) - spin_field(get_south_neighbour(i)) &
                        + n * twopi / integer_lattice_length + pi, twopi) - pi) ** 2
    end do
    potential_difference = 0.5d0 * (twisted_sum_of_squared_electric_field_x - current_sum_of_squared_electric_field_x)
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential_y = n
        n = n + 1
    else
        exit
    end if
end do

! y direction; negative twist
n = 1
do
    current_sum_of_squared_electric_field_x = sum_of_squared_electric_field_x
    twisted_sum_of_squared_electric_field_x = 0.0d0
    do i = 1, no_of_sites
        twisted_sum_of_squared_electric_field_x = twisted_sum_of_squared_electric_field_x &
                + (modulo(spin_field(i) - spin_field(get_south_neighbour(i)) &
                        - n * twopi / integer_lattice_length + pi, twopi) - pi) ** 2
    end do
    potential_difference = 0.5d0 * (twisted_sum_of_squared_electric_field_x - current_sum_of_squared_electric_field_x)
    if (potential_difference < 0.0d0) then
        no_of_external_twists_to_minimise_potential_y = n
        n = n + 1
    else
        exit
    end if
end do

return
end subroutine
