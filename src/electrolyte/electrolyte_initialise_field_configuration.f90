subroutine initialise_field_configuration
use variables
implicit none
integer :: i

do i = 1, no_of_sites
    electric_field_x(i) = 0.0
    electric_field_y(i) = 0.0
end do

electric_field_sum_x = 0.0d0
electric_field_sum_y = 0.0d0

return
end subroutine initialise_field_configuration
