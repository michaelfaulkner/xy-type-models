! **************************************
! MEASURE STATE OF SYSTEM
! **************************************

subroutine measure
  use variables
  implicit none
  integer i, n
  real*8 magn, magn_x, magn_y, cos_top_x, cos_top_y
  real*8 sin_top_x, sin_top_y, Ebar_x, Ebar_y, vort
  real*8 potential, storeTop_x, storeTop_y

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
  potential = 0.0

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
     potential = potential - cos(storeTop_x) - cos(storeTop_y)

  end do


  magn = dsqrt(magn_x ** 2 + magn_y ** 2)
  magn = magn / volume
  magn_x = magn_x / volume
  magn_y = magn_y / volume
  cos_top_x = cos_top_x / volume
  cos_top_y = cos_top_y / volume
  sin_top_x = sin_top_x / volume
  sin_top_y = sin_top_y / volume
  Ebar_x = Ebar_x / volume
  Ebar_y = Ebar_y / volume
  vort = vort / volume
  

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
  write(21, 100) potential

100 format(F16.8)

  return
end subroutine measure


! **************************************
! CALCULATE EMERGENT FIELD
! VIA EMERGENT-FIELD DEF: E_x(i) = + (theta(i) - theta(neg_y(i))) (WITH MODULAR OPERATION; SIMILARLY IN Y, UP TO A MINUS SIGN)
! **************************************

subroutine top_field
use variables
implicit none
real*8 diff
integer i, j

do i = 1, sites
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
  do i = 1, sites
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
