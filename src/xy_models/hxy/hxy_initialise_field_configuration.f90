subroutine initialise_field_configuration(pre_simulation)
use variables
implicit none
logical :: pre_simulation
integer :: i

if (randomise_initial_field_configuration) then
    do i = 1, no_of_sites
        spin_field(i) = two_pi * rand()
    end do
    call calculate_emergent_field
    external_global_move = (/ 0, 0 /)
else if (pre_simulation) then
    do i = 1, no_of_sites
        spin_field(i) = 0.0d0
    end do
    call calculate_emergent_field
    external_global_move = (/ 0, 0 /)
end if

return
end subroutine initialise_field_configuration
