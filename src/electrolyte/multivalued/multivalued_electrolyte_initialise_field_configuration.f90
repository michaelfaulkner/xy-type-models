subroutine initialise_field_configuration(pre_simulation)
use variables
implicit none
logical :: pre_simulation
integer :: i

if (pre_simulation) then
    do i = 1, no_of_sites
        electric_field(i, 1) = 0.0
        electric_field(i, 2) = 0.0
    end do
    net_charge_displacement = (/ 0, 0 /)
end if

return
end subroutine initialise_field_configuration
