subroutine draw_observations
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

do i = 1, sites
    magnetisation_x = magnetisation_x + cos(theta(i))
    magnetisation_y = magnetisation_y + sin(theta(i))
    sum_of_1st_derivative_of_potential_x = sum_of_1st_derivative_of_potential_x + sin(theta(pos_x(i)) - theta(i))
    sum_of_1st_derivative_of_potential_y = sum_of_1st_derivative_of_potential_y + sin(theta(pos_y(i)) - theta(i))
    sum_of_2nd_derivative_of_potential_x = sum_of_2nd_derivative_of_potential_x + cos(theta(pos_x(i)) - theta(i))
    sum_of_2nd_derivative_of_potential_y = sum_of_2nd_derivative_of_potential_y + cos(theta(pos_y(i)) - theta(i))
    potential = potential - cos(theta(pos_x(i)) - theta(i)) - cos(theta(pos_y(i)) - theta(i))
end do

write(20, 100) potential, magnetisation_x, magnetisation_y, &
                sum_of_1st_derivative_of_potential_x, sum_of_1st_derivative_of_potential_y, &
                sum_of_2nd_derivative_of_potential_x, sum_of_2nd_derivative_of_potential_y

100 format(ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14, ",", ES29.14)

return
end subroutine draw_observations
