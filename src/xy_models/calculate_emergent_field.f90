subroutine calculate_emergent_field
use variables
implicit none
integer :: i
double precision :: get_spin_difference

do i = 1, no_of_sites
    emergent_field(i, 1) = get_spin_difference(spin_field(i), spin_field(get_south_neighbour(i)))
    emergent_field(i, 2) = get_spin_difference(spin_field(get_west_neighbour(i)), spin_field(i))
end do

return
end subroutine calculate_emergent_field