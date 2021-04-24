subroutine calculate_emergent_field
use variables
implicit none
integer :: i
do i = 1, no_of_sites
    emergent_field(i, 1) = modulo(spin_field(i) - spin_field(get_south_neighbour(i)) + pi, twopi) - pi
    emergent_field(i, 2) = modulo(- spin_field(i) + spin_field(get_west_neighbour(i)) + pi, twopi) - pi
end do
return
end subroutine calculate_emergent_field