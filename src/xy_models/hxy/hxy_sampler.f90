! **************************************
! MEASURE STATE OF SYSTEM
! **************************************

subroutine draw_observations
use variables
implicit none
integer :: i, n
double precision :: magnetisation_x, magnetisation_y, potential
double precision :: mean_1st_derivative_of_potential_x, mean_1st_derivative_of_potential_y
double precision :: mean_2nd_derivative_of_potential_x, mean_2nd_derivative_of_potential_y

call top_field

magnetisation_x = 0.0d0
magnetisation_y = 0.0d0
mean_1st_derivative_of_potential_x = 0.0d0
mean_1st_derivative_of_potential_y = 0.0d0
mean_2nd_derivative_of_potential_x = 0.0d0
mean_2nd_derivative_of_potential_y = 0.0d0
sum_of_squared_electric_field_x = 0.0d0
sum_of_squared_electric_field_y = 0.0d0

do i = 1, sites
    magnetisation_x = magnetisation_x + cos(theta(i))
    magnetisation_y = magnetisation_y + sin(theta(i))
    mean_1st_derivative_of_potential_x = mean_1st_derivative_of_potential_x + top_y(i)
    mean_1st_derivative_of_potential_y = mean_1st_derivative_of_potential_y + top_x(i)
    do n = 1, nmax
        mean_2nd_derivative_of_potential_x = mean_2nd_derivative_of_potential_x + (-1) ** (n + 1) * cos((n + 1) * top_y(i))
        mean_2nd_derivative_of_potential_y = mean_2nd_derivative_of_potential_y + (-1) ** (n + 1) * cos((n + 1) * top_x(i))
    end do
    sum_of_squared_electric_field_x = sum_of_squared_electric_field_x + top_x(i) * top_x(i)
    sum_of_squared_electric_field_y = sum_of_squared_electric_field_y + top_y(i) * top_y(i)
end do
potential = 0.5d0 * (sum_of_squared_electric_field_x + sum_of_squared_electric_field_y)

magnetisation_x = magnetisation_x / volume
magnetisation_y = magnetisation_y / volume
mean_1st_derivative_of_potential_x = mean_1st_derivative_of_potential_x / volume
mean_1st_derivative_of_potential_y = mean_1st_derivative_of_potential_y / volume
mean_2nd_derivative_of_potential_x = mean_2nd_derivative_of_potential_x / volume
mean_2nd_derivative_of_potential_y = mean_2nd_derivative_of_potential_y / volume

if (calculate_external_minimising_twist_field == 1) then
    call external_minimising_twist_field_calculation
end if
  
write(10, 100) magnetisation_x
write(11, 100) magnetisation_y
write(12, 100) mean_1st_derivative_of_potential_x
write(13, 100) mean_1st_derivative_of_potential_y
write(14, 100) mean_2nd_derivative_of_potential_x
write(15, 100) mean_2nd_derivative_of_potential_y
write(16, 100) potential
write(17,200) no_of_external_twists_to_minimise_potential_x
write(18,200) no_of_external_twists_to_minimise_potential_y

100 format(ES24.14)
200 format(I2)

return
end subroutine draw_observations

! **************************************
! CALCULATE EMERGENT FIELD
! VIA EMERGENT-FIELD DEF: E_x(i) = + (theta(i) - theta(neg_y(i)));
! E_y(i) = - (theta(i) - theta(neg_x(i))) (WITH MODULAR OPERATION)
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


! **************************************
! Calculate the number of external twists required to minimise the current value of the potential
! **************************************

subroutine external_minimising_twist_field_calculation
    use variables
    implicit none
    integer i, n
    real*8 diff, current_sum_of_squared_electric_field_x, current_sum_of_squared_electric_field_y
    real*8 twisted_sum_of_squared_electric_field_x, twisted_sum_of_squared_electric_field_y, potential_difference

  ! y direction; positive twist
  n = 1
  do
     current_sum_of_squared_electric_field_x = sum_of_squared_electric_field_x
     twisted_sum_of_squared_electric_field_x = 0.0

     do i = 1, sites
        diff = theta(i) - theta(neg_y(i)) + n * twopi / side
        if (diff .gt. 0.5 * twopi) then
           diff = diff - twopi
        else if (diff .le. - 0.5 * twopi) then
           diff = diff + twopi
        end if
        twisted_sum_of_squared_electric_field_x = twisted_sum_of_squared_electric_field_x + diff ** 2
     end do

     potential_difference = 0.5 * (twisted_sum_of_squared_electric_field_x - current_sum_of_squared_electric_field_x)

     if (potential_difference .lt. - epsilon) then
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
     twisted_sum_of_squared_electric_field_x = 0.0

     do i = 1, sites
        diff = theta(i) - theta(neg_y(i)) - n * twopi / side
        if (diff .gt. 0.5 * twopi) then
           diff = diff - twopi
        else if (diff .le. - 0.5 * twopi) then
           diff = diff + twopi
        end if
        twisted_sum_of_squared_electric_field_x = twisted_sum_of_squared_electric_field_x + diff ** 2
     end do

     potential_difference = 0.5 * (twisted_sum_of_squared_electric_field_x - current_sum_of_squared_electric_field_x)

     if (potential_difference .lt. - epsilon) then
         no_of_external_twists_to_minimise_potential_y = - n
        n = n + 1
     else
        exit
     end if
  end do

  ! x direction; positive twist
  n = 1
  do
     current_sum_of_squared_electric_field_y = sum_of_squared_electric_field_y
     twisted_sum_of_squared_electric_field_y = 0.0

     do i = 1, sites
        diff = - (theta(i) - theta(neg_x(i))) + n * twopi / side
        if (diff .gt. 0.5 * twopi) then
           diff = diff - twopi
        else if (diff .le. - 0.5 * twopi) then
           diff = diff + twopi
        end if
        twisted_sum_of_squared_electric_field_y = twisted_sum_of_squared_electric_field_y + diff ** 2
     end do

     potential_difference = 0.5 * (twisted_sum_of_squared_electric_field_y - current_sum_of_squared_electric_field_y)

     if (potential_difference .lt. - epsilon) then
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
     twisted_sum_of_squared_electric_field_y = 0.0

     do i = 1, sites
        diff = - (theta(i) - theta(neg_x(i))) - n * twopi / side
        if (diff .gt. 0.5 * twopi) then
           diff = diff - twopi
        else if (diff .le. - 0.5 * twopi) then
           diff = diff + twopi
        end if
        twisted_sum_of_squared_electric_field_y = twisted_sum_of_squared_electric_field_y + diff ** 2
     end do

     potential_difference = 0.5 * (twisted_sum_of_squared_electric_field_y - current_sum_of_squared_electric_field_y)

     if (potential_difference .lt. - epsilon) then
         no_of_external_twists_to_minimise_potential_x = - n
        n = n + 1
     else
        exit
     end if
  end do

end subroutine
