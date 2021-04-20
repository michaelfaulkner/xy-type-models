subroutine initialise_field_configuration
use variables
implicit none
integer :: i

do i = 1, sites
    Efield_x(i) = 0.0
    Efield_y(i) = 0.0
end do

return
end subroutine initialise_field_configuration
