subroutine calculate_charge_configuration
use variables
implicit none
integer :: i

do i = 1, sites
    rho(i) = floor((electric_field_x(i) + electric_field_y(i) - electric_field_x(neg_x(i)) - &
                        electric_field_y(neg_y(i))) / twopi)
end do

return
end subroutine calculate_charge_configuration
