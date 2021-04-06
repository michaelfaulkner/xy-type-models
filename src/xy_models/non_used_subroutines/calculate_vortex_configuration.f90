subroutine calculate_vortex_configuration
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
end subroutine calculate_vortex_configuration