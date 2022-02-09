subroutine get_and_print_observation
use variables
implicit none
integer :: i
double precision :: potential, non_normalised_magnetisation(2)
double precision :: sum_of_1st_derivative_of_potential(2), sum_of_2nd_derivative_of_potential(2)

potential = 0.0d0
non_normalised_magnetisation = (/ 0.0d0, 0.0d0 /)
sum_of_1st_derivative_of_potential = (/ 0.0d0, 0.0d0 /)
sum_of_2nd_derivative_of_potential = (/ 0.0d0, 0.0d0 /)

do i = 1, no_of_sites
    non_normalised_magnetisation(1) = non_normalised_magnetisation(1) + cos(spin_field(i))
    non_normalised_magnetisation(2) = non_normalised_magnetisation(2) + sin(spin_field(i))
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

write(20, 100) potential, external_global_move(1), external_global_move(2), &
                non_normalised_magnetisation(1), non_normalised_magnetisation(2), &
                sum_of_1st_derivative_of_potential(1), sum_of_1st_derivative_of_potential(2), &
                sum_of_2nd_derivative_of_potential(1), sum_of_2nd_derivative_of_potential(2)

100 format(ES30.14, ",", I29, ",", I29, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", &
            ES29.14)

return
end subroutine get_and_print_observation
