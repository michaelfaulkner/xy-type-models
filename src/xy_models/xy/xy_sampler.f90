subroutine draw_observations
use variables
implicit none
integer :: i
double precision :: magnetisation_x, magnetisation_y, potential
double precision :: mean_1st_derivative_of_potential_x, mean_1st_derivative_of_potential_y
double precision :: mean_2nd_derivative_of_potential_x, mean_2nd_derivative_of_potential_y

magnetisation_x = 0.0d0
magnetisation_y = 0.0d0
mean_1st_derivative_of_potential_x = 0.0d0
mean_1st_derivative_of_potential_y = 0.0d0
mean_2nd_derivative_of_potential_x = 0.0d0
mean_2nd_derivative_of_potential_y = 0.0d0
potential = 0.0d0

do i = 1, sites
    magnetisation_x = magnetisation_x + cos(theta(i))
    magnetisation_y = magnetisation_y + sin(theta(i))
    mean_1st_derivative_of_potential_x = mean_1st_derivative_of_potential_x + sin(theta(pos_x(i)) - theta(i))
    mean_1st_derivative_of_potential_y = mean_1st_derivative_of_potential_y + sin(theta(pos_y(i)) - theta(i))
    mean_2nd_derivative_of_potential_x = mean_2nd_derivative_of_potential_x + cos(theta(pos_x(i)) - theta(i))
    mean_2nd_derivative_of_potential_y = mean_2nd_derivative_of_potential_y + cos(theta(pos_y(i)) - theta(i))
    potential = potential - cos(theta(pos_x(i)) - theta(i)) - cos(theta(pos_y(i)) - theta(i))
end do

magnetisation_x = magnetisation_x / volume
magnetisation_y = magnetisation_y / volume
mean_1st_derivative_of_potential_x = mean_1st_derivative_of_potential_x / volume
mean_1st_derivative_of_potential_y = mean_1st_derivative_of_potential_y / volume
mean_2nd_derivative_of_potential_x = mean_2nd_derivative_of_potential_x / volume
mean_2nd_derivative_of_potential_y = mean_2nd_derivative_of_potential_y / volume

write(10, 100) magnetisation_x, magnetisation_y
write(11, 100) mean_1st_derivative_of_potential_x, mean_1st_derivative_of_potential_y
write(12, 100) mean_2nd_derivative_of_potential_x, mean_2nd_derivative_of_potential_y
write(13, 200) potential

100 format(ES24.14, ", ", ES24.14)
200 format(ES24.14)

return
end subroutine draw_observations
