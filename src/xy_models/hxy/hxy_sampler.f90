subroutine draw_observations
use variables
implicit none
integer :: i, n
double precision :: magnetisation_x, magnetisation_y, potential
double precision :: sum_of_1st_derivative_of_potential_x, sum_of_1st_derivative_of_potential_y
double precision :: sum_of_2nd_derivative_of_potential_x, sum_of_2nd_derivative_of_potential_y

call top_field

magnetisation_x = 0.0d0
magnetisation_y = 0.0d0
sum_of_1st_derivative_of_potential_x = 0.0d0
sum_of_1st_derivative_of_potential_y = 0.0d0
sum_of_2nd_derivative_of_potential_x = 0.0d0
sum_of_2nd_derivative_of_potential_y = 0.0d0
sum_of_squared_electric_field_x = 0.0d0
sum_of_squared_electric_field_y = 0.0d0

do i = 1, sites
    magnetisation_x = magnetisation_x + cos(theta(i))
    magnetisation_y = magnetisation_y + sin(theta(i))
    sum_of_1st_derivative_of_potential_x = sum_of_1st_derivative_of_potential_x + top_y(i)
    sum_of_1st_derivative_of_potential_y = sum_of_1st_derivative_of_potential_y + top_x(i)
    do n = 1, nmax
        sum_of_2nd_derivative_of_potential_x = sum_of_2nd_derivative_of_potential_x + &
                                                                                (-1.0d0) ** (n + 1) * cos(n * top_y(i))
        sum_of_2nd_derivative_of_potential_y = sum_of_2nd_derivative_of_potential_y + &
                                                                                (-1.0d0) ** (n + 1) * cos(n * top_x(i))
    end do
    sum_of_squared_electric_field_x = sum_of_squared_electric_field_x + top_x(i) * top_x(i)
    sum_of_squared_electric_field_y = sum_of_squared_electric_field_y + top_y(i) * top_y(i)
end do
sum_of_2nd_derivative_of_potential_x = 2.0d0 * sum_of_2nd_derivative_of_potential_x
sum_of_2nd_derivative_of_potential_y = 2.0d0 * sum_of_2nd_derivative_of_potential_y
potential = 0.5d0 * (sum_of_squared_electric_field_x + sum_of_squared_electric_field_y)

if (calculate_external_minimising_twist_field) then
    call external_minimising_twist_field_calculation
end if
  
write(20, 100) potential, magnetisation_x, magnetisation_y, &
                sum_of_1st_derivative_of_potential_x, sum_of_1st_derivative_of_potential_y, &
                sum_of_2nd_derivative_of_potential_x, sum_of_2nd_derivative_of_potential_y, &
                no_of_external_twists_to_minimise_potential_x, no_of_external_twists_to_minimise_potential_y

100 format(ES30.14, ", ", ES30.14, ", ", ES30.14, ", ", ES30.14, ", ", ES30.14, ", ", ES30.14, ", ", ES30.14, ", ", &
            I2, ", ", I2)

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
! Calculate the number of external twists required to minimise the current value of the potential
! **************************************

subroutine external_minimising_twist_field_calculation
    use variables
    implicit none
    integer i, n
    double precision diff, current_sum_of_squared_electric_field_x, current_sum_of_squared_electric_field_y
    double precision twisted_sum_of_squared_electric_field_x, twisted_sum_of_squared_electric_field_y, potential_difference

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
