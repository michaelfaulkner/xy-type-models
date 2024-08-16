subroutine initialise_field_configuration(pre_simulation)
use variables
implicit none
logical :: pre_simulation
integer :: site_index

if ((pre_simulation).or.(always_cold_start)) then
    do site_index = 1, no_of_sites
        electric_field(site_index, 1) = 0.0d0
        electric_field(site_index, 2) = 0.0d0
        charge_configuration(site_index) = 0
    end do
    net_charge_displacement = (/ 0, 0 /)
    external_global_moves = (/ 0, 0 /)
end if

return
end subroutine initialise_field_configuration
