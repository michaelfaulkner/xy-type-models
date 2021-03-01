! **************************************
! MEASURE STATE OF SYSTEM
! **************************************

subroutine measure
  use variables
  implicit none
  integer i, j, n
  real*8 magn, magn_x, magn_y, cos_top_x, cos_top_y, sin_top_x, sin_top_y
  real*8 potential, vort, cos_j, second_deriv_potential, diff, untwist, Ebar_x, Ebar_y
  real*8 storeTop_x, storeTop_y

  call top_field
  call vortices

  magn_x = 0.0
  magn_y = 0.0
  cos_top_x = 0.0
  cos_top_y = 0.0
  sin_top_x = 0.0
  sin_top_y = 0.0
  Ebar_x = 0.0
  Ebar_y = 0.0
  vort = 0.0
  second_deriv_potential = 0.0
  potential = 0.0

  do i = 0, sites - 1
     
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
     potential = potential + 0.5 * (storeTop_x * storeTop_x + storeTop_y * storeTop_y)
     
  end do

  ! CALCULATE THE CUT-OFF FOURIER SERIES THAT APPROXIMATES THE SECOND DERIVATIVE OF THE POTENTIAL
  ! THIS IS THE FIRST TERM OF THE HELICITY MODULUS, AND IS PRECISE FOR nmax = infty
  
  do j = 0, nmax - 1
     cos_j = 0.0
     do i = 0, sites - 1
        cos_j = cos_j + cos((j + 1) * top_x(i)) + cos((j + 1) * top_y(i))
     end do
     second_deriv_potential = second_deriv_potential + (-1) ** (j + 2) * cos_j
  end do

  magn = dsqrt(magn_x ** 2 + magn_y ** 2)
  magn = magn / volume
  magn_x = magn_x / volume
  magn_y = magn_y / volume
  cos_top_x = cos_top_x / volume
  cos_top_y = cos_top_y / volume
  sin_top_x = sin_top_x / volume
  sin_top_y = sin_top_y / volume
  vort = vort / volume
  second_deriv_potential = second_deriv_potential / volume
  
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

100 format(F16.8)

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
  real*8 diff
  integer i,j
  do i = 0, sites - 1

     diff = theta(i) - theta(neg_y(i))
     if (diff .gt. 0.5 * twopi) then
        diff = diff - twopi
     else if (diff .le. -0.5 * twopi) then
        diff = diff + twopi
     end if
     top_x(i) = diff
     
     diff = -(theta(i) - theta(neg_x(i)))
     if (diff .gt. 0.5 * twopi) then
        diff = diff - twopi
     else if (diff .le. -0.5 * twopi) then
        diff = diff + twopi
     end if
     top_y(i) = diff
     
  end do
  return
end subroutine top_field


! **************************************
! MEASURE NUMBER OF VORTICES
! **************************************

subroutine vortices
  use variables
  implicit none
  integer i,j
  real*8 v0
  do i = 0, sites - 1
     v0 = (top_x(i) + top_y(i) - top_x(neg_x(i)) - top_y(neg_y(i))) / twopi
     if ((v0 .lt. 1.001) .and. (v0 .gt. 0.999)) then
        v(i) = 1
     else if ((v0 .gt. -1.001) .and. (v0 .lt. -0.999)) then
        v(i) = -1
     else
        v(i) = 0
     end if
  end do
  return
end subroutine vortices


! **************************************
! OUTPUT acceptance rates
! **************************************

subroutine output_acceptance_rates
  use variables
  implicit none
  character(100) filename
  real*8 acceptanceRate, twistAcceptanceRate

  acceptanceRate = float(accept) / (measurements * volume)
  twistAcceptanceRate = float(accept_twist) / (measurements * volume)

  write (filename, '( "temp_eq_", F4.2,"//acceptance_rates_HXY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' ) &
       T, side, side, T
  open(unit = 300, file = filename)
  
  write(300, 100) acceptanceRate
  if (twist .eq. 1) then
     write(300, 100) twistAcceptanceRate
  end if
  close(300)

100 format(F16.8)
  
  return
end subroutine output_acceptance_rates
