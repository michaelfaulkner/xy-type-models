subroutine calculate_charge_configuration
use variables
implicit none
integer :: i

do i = 1, no_of_sites
    rho(i) = floor((electric_field_x(i) + electric_field_y(i) - electric_field_x(get_west_neighbour(i)) - &
                        electric_field_y(get_south_neighbour(i))) / twopi)
end do

return
end subroutine calculate_charge_configuration
