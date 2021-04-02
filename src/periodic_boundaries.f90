subroutine setup_periodic_boundaries
use variables
implicit none
integer :: i

! mod(i + side + 1, side) rather than mod(i + 1, side) below as mod(i, side) doesn't return value in [0, side) if i < 0
do i = 1, sites
    pos_x(i) = i + mod(i + side + 1, side) - mod(i + side, side)
    neg_x(i) = i + mod(i + side - 1, side) - mod(i + side, side)
    pos_y(i) = i + (mod(int(i / side) + side + 1, side) - mod(int(i / side) + side, side)) * side
    neg_y(i) = i + (mod(int(i / side) + side - 1, side) - mod(int(i / side) + side, side)) * side
end do

return
end subroutine setup_periodic_boundaries
