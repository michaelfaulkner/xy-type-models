subroutine calculate_emergent_field
use variables
implicit none
integer :: i
do i = 1, no_of_sites
    emergent_field_x(i) = modulo(spin_field(i) - spin_field(neg_y(i)) + pi, twopi) - pi
    emergent_field_y(i) = modulo(- spin_field(i) + spin_field(neg_x(i)) + pi, twopi) - pi
end do
return
end subroutine calculate_emergent_field