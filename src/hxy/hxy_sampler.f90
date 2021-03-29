! **************************************
! MEASURE STATE OF SYSTEM
! **************************************

subroutine measure
  use variables
  implicit none
  integer :: i, j, n
  double precision :: magn, magn_x, magn_y, cos_top_x, cos_top_y, sin_top_x, sin_top_y
  double precision :: potential, vort, cos_j, second_deriv_potential, diff, Ebar_x, Ebar_y, storeTop_x, storeTop_y

  call top_field
  call vortices

  magn_x = 0.0d0
  magn_y = 0.0d0
  cos_top_x = 0.0d0
  cos_top_y = 0.0d0
  sin_top_x = 0.0d0
  sin_top_y = 0.0d0
  Ebar_x = 0.0d0
  Ebar_y = 0.0d0
  vort = 0.0d0
  second_deriv_potential = 0.0d0

  sum_of_squared_electric_field_x = 0.0d0
  sum_of_squared_electric_field_y = 0.0d0

  do i = 1, sites
     
     magn_x = magn_x + cos(theta(i))
     magn_y = magn_y + sin(theta(i))
     storeTop_x = top_x(i)
     storeTop_y = top_y(i)
     cos_top_x = cos_top_x + cos(storeTop_x)
     cos_top_y = cos_top_y + cos(storeTop_y)
     sin_top_x = sin_top_x + sin(storeTop_x)
     sin_top_y = sin_top_y + sin(storeTop_y)
     Ebar_x = Ebar_x + storeTop_x
     Ebar_y = Ebar_y + storeTop_y
     if (v(i) .ne. 0) then
        vort = vort + 1.0
     end if

     sum_of_squared_electric_field_x = sum_of_squared_electric_field_x + storeTop_x ** 2
     sum_of_squared_electric_field_y = sum_of_squared_electric_field_y + storeTop_y ** 2
     
  end do

  potential = 0.5 * (sum_of_squared_electric_field_x + sum_of_squared_electric_field_y)
  
  ! CALCULATE THE CUT-OFF FOURIER SERIES THAT APPROXIMATES THE SECOND DERIVATIVE OF THE POTENTIAL
  ! THIS IS THE FIRST TERM OF THE HELICITY MODULUS, AND IS PRECISE FOR nmax = infty

  do j = 1, nmax
     cos_j = 0.0d0
     do i = 1, sites
        cos_j = cos_j + cos((j + 1) * top_x(i)) + cos((j + 1) * top_y(i))
     end do
     second_deriv_potential = second_deriv_potential + (-1) ** (j + 2) * cos_j
  end do

  magn = sqrt(magn_x ** 2 + magn_y ** 2)
  magn = magn / volume
  magn_x = magn_x / volume
  magn_y = magn_y / volume
  cos_top_x = cos_top_x / volume
  cos_top_y = cos_top_y / volume
  sin_top_x = sin_top_x / volume
  sin_top_y = sin_top_y / volume
  vort = vort / volume
  second_deriv_potential = second_deriv_potential / volume
  Ebar_x = Ebar_x / volume
  Ebar_y = Ebar_y / volume

    if (calculate_external_minimising_twist_field == 1) then
        call external_minimising_twist_field_calculation
    end if
  
  ! STORE SAMPLES DRAWN FROM MARKOV CHAIN
  
  write(10,100) magn
  write(11,100) magn_x
  write(12,100) magn_y
  write(13,100) cos_top_x
  write(14,100) cos_top_y
  write(15,100) sin_top_x
  write(16,100) sin_top_y
  write(17,100) Ebar_x
  write(18,100) Ebar_y
  write(19,100) vort
  write(20,100) second_deriv_potential
  write(21,100) potential

  write(22,200) no_of_external_twists_to_minimise_potential_x
  write(23,200) no_of_external_twists_to_minimise_potential_y


100 format(ES24.14)
200 format(I2)

  return
end subroutine measure

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
