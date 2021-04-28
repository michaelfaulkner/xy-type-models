subroutine initialise_field_configuration
use variables
implicit none
integer :: i

do i = 1, no_of_sites
    electric_field(i, 1) = 0.0
    electric_field(i, 2) = 0.0
    charge_configuration(i) = 0
end do
net_charge_displacement = (/ 0.0d0, 0.0d0 /)

return
end subroutine initialise_field_configuration
