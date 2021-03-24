! **************************************
! MEASURE STATE OF SYSTEM
! **************************************

subroutine measure
use variables
implicit none
integer :: i, n
double precision :: magn, magn_x, magn_y, cos_top_x, cos_top_y, sin_top_x, sin_top_y
double precision :: Ebar_x, Ebar_y, vort, potential, storeTop_x, storeTop_y

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
potential = 0.0d0

do i = 1, sites
    magn_x = magn_x + dcos(theta(i))
    magn_y = magn_y + dsin(theta(i))
    storeTop_x = top_x(i)
    storeTop_y = top_y(i)
    cos_top_x = cos_top_x + dcos(storeTop_x)
    cos_top_y = cos_top_y + dcos(storeTop_y)
    sin_top_x = sin_top_x + dsin(storeTop_x)
    sin_top_y = sin_top_y + dsin(storeTop_y)
    Ebar_x = Ebar_x + storeTop_x
    Ebar_y = Ebar_y + storeTop_y
    if (v(i) /= 0) then
        vort = vort + 1.0
    end if
    potential = potential - dcos(storeTop_x) - dcos(storeTop_y)
end do

magn = dsqrt(magn_x ** 2 + magn_y ** 2) / volume
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

write(10, 100) magn
write(11, 100) magn_x
write(12, 100) magn_y
write(13, 100) cos_top_x
write(14, 100) cos_top_y
write(15, 100) sin_top_x
write(16, 100) sin_top_y
write(17, 100) Ebar_x
write(18, 100) Ebar_y
write(19, 100) vort
write(21, 100) potential

100 format(ES24.14)

return
end subroutine measure


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
    if ((v0 < 1.000000000001) .and. (v0 > 0.999999999999)) then
        v(i) = 1
    else if ((v0 > -1.000000000001) .and. (v0 < -0.999999999999)) then
        v(i) = -1
    else
        v(i) = 0
    end if
end do
return
end subroutine vortices
