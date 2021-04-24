subroutine get_and_print_observation
use variables
implicit none
integer :: i
double precision :: potential
double precision :: magnetisation(2), sum_of_1st_derivative_of_potential(2), sum_of_2nd_derivative_of_potential(2)

potential = 0.0d0
magnetisation = (/ 0.0d0, 0.0d0 /)
sum_of_1st_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
sum_of_2nd_derivative_of_potential = (/ 0.0d0, 0.0d0 /)

do i = 1, no_of_sites
    magnetisation(1) = magnetisation(1) + cos(spin_field(i))
    magnetisation(2) = magnetisation(2) + sin(spin_field(i))
    sum_of_1st_derivative_of_potential(1) = sum_of_1st_derivative_of_potential(1) &
                                            + sin(spin_field(get_east_neighbour(i)) - spin_field(i))
    sum_of_1st_derivative_of_potential(2) = sum_of_1st_derivative_of_potential(2) &
                                            + sin(spin_field(get_north_neighbour(i)) - spin_field(i))
    sum_of_2nd_derivative_of_potential(1) = sum_of_2nd_derivative_of_potential(1) &
                                            + cos(spin_field(get_east_neighbour(i)) - spin_field(i))
    sum_of_2nd_derivative_of_potential(2) = sum_of_2nd_derivative_of_potential(2) &
                                            + cos(spin_field(get_north_neighbour(i)) - spin_field(i))
    potential = potential - cos(spin_field(get_east_neighbour(i)) - spin_field(i)) &
                          - cos(spin_field(get_north_neighbour(i)) - spin_field(i))
end do

write(20, 100) potential, magnetisation(1), magnetisation(2), &
                sum_of_1st_derivative_of_potential(1), sum_of_1st_derivative_of_potential(2), &
                sum_of_2nd_derivative_of_potential(1), sum_of_2nd_derivative_of_potential(2)

100 format(ES30.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14)

return
end subroutine get_and_print_observation
