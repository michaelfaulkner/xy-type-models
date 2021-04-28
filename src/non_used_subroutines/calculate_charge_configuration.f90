subroutine calculate_charge_configuration
use variables
implicit none
integer :: i

do i = 1, no_of_sites
    rho(i) = floor((electric_field(i, 1) + electric_field(i, 2)(i) - electric_field(get_west_neighbour(i), 1) - &
                        electric_field(get_south_neighbour(i), 2)) / two_pi)
end do

return
end subroutine calculate_charge_configuration
