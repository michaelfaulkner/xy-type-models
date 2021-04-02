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

write(10, 100) magnetisation_x
write(11, 100) magnetisation_y
write(12, 100) mean_1st_derivative_of_potential_x
write(13, 100) mean_1st_derivative_of_potential_y
write(14, 100) mean_2nd_derivative_of_potential_x
write(15, 100) mean_2nd_derivative_of_potential_y
write(16, 100) potential

100 format(ES24.14)

return
end subroutine draw_observations


! **************************************
! CALCULATE EMERGENT FIELD
! VIA EMERGENT-FIELD DEF: E_x(i) = + (theta(i) - theta(neg_y(i))) (WITH MODULAR OPERATION; SIMILARLY IN Y, UP TO A MINUS SIGN)
! **************************************

subroutine top_field
use variables
implicit none
integer :: i
do i = 1, sites
    top_x(i) = modulo(theta(i) - theta(neg_y(i)) + pi, twopi) - pi
    top_y(i) = modulo(- theta(i) + theta(neg_x(i)) + pi, twopi) - pi
end do
return
end subroutine top_field


! **************************************
! MEASURE NUMBER OF VORTICES
! **************************************

subroutine vortices
use variables
implicit none
integer :: i
double precision :: v0
do i = 1, sites
    v0 = (top_x(i) + top_y(i) - top_x(neg_x(i)) - top_y(neg_y(i))) / twopi
    if ((v0 < 1.000000000001d0) .and. (v0 > 0.999999999999d0)) then
        v(i) = 1
    else if ((v0 > -1.000000000001d0) .and. (v0 < -0.999999999999d0)) then
        v(i) = -1
    else
        v(i) = 0
    end if
end do
return
end subroutine vortices
