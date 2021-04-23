subroutine get_and_print_observation
use variables
implicit none
integer :: i
double precision :: magnetisation_x, magnetisation_y, potential
double precision :: sum_of_1st_derivative_of_potential_x, sum_of_1st_derivative_of_potential_y
double precision :: sum_of_2nd_derivative_of_potential_x, sum_of_2nd_derivative_of_potential_y

magnetisation_x = 0.0d0
magnetisation_y = 0.0d0
sum_of_1st_derivative_of_potential_x = 0.0d0
sum_of_1st_derivative_of_potential_y = 0.0d0
sum_of_2nd_derivative_of_potential_x = 0.0d0
sum_of_2nd_derivative_of_potential_y = 0.0d0
potential = 0.0d0

do i = 1, no_of_sites
    magnetisation_x = magnetisation_x + cos(spin_field(i))
    magnetisation_y = magnetisation_y + sin(spin_field(i))
    sum_of_1st_derivative_of_potential_x = sum_of_1st_derivative_of_potential_x &
                                            + sin(spin_field(get_east_neighbour(i)) - spin_field(i))
    sum_of_1st_derivative_of_potential_y = sum_of_1st_derivative_of_potential_y &
                                            + sin(spin_field(get_north_neighbour(i)) - spin_field(i))
    sum_of_2nd_derivative_of_potential_x = sum_of_2nd_derivative_of_potential_x &
                                            + cos(spin_field(get_east_neighbour(i)) - spin_field(i))
    sum_of_2nd_derivative_of_potential_y = sum_of_2nd_derivative_of_potential_y &
                                            + cos(spin_field(get_north_neighbour(i)) - spin_field(i))
    potential = potential - cos(spin_field(get_east_neighbour(i)) - spin_field(i)) &
                          - cos(spin_field(get_north_neighbour(i)) - spin_field(i))
end do

write(20, 100) potential, magnetisation_x, magnetisation_y, &
                sum_of_1st_derivative_of_potential_x, sum_of_1st_derivative_of_potential_y, &
                sum_of_2nd_derivative_of_potential_x, sum_of_2nd_derivative_of_potential_y

100 format(ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14)

return
end subroutine get_and_print_observation
