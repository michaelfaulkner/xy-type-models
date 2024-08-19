subroutine initialise_field_configuration(pre_simulation)
use variables
implicit none
logical :: pre_simulation
integer :: site_index

if (always_hot_start) then
    do site_index = 1, no_of_sites
        spin_field(site_index) = two_pi * rand()
    end do
    external_global_moves = (/ 0, 0 /)
else if ((always_cold_start).or.(pre_simulation)) then
    do site_index = 1, no_of_sites
        spin_field(site_index) = 0.0d0
    end do
    external_global_moves = (/ 0, 0 /)
end if

return
end subroutine initialise_field_configuration
