subroutine calculate_emergent_field
use variables
implicit none
integer :: i
do i = 1, sites
    emergent_field_x(i) = modulo(theta(i) - theta(neg_y(i)) + pi, twopi) - pi
    emergent_field_y(i) = modulo(- theta(i) + theta(neg_x(i)) + pi, twopi) - pi
end do
return
end subroutine calculate_emergent_field