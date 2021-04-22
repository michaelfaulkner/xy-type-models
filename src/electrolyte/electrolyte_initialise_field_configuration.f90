subroutine initialise_field_configuration
use variables
implicit none
integer :: i

do i = 1, sites
    electric_field_x(i) = 0.0
    electric_field_y(i) = 0.0
end do

return
end subroutine initialise_field_configuration
