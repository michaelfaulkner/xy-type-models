subroutine initialise_field_configuration
use variables
implicit none
integer :: i

if (start == 0) then
    do i = 1, sites
        theta(i) = 0.0d0
    end do
else
    do i = 1, sites
        theta(i) = twopi * rand()
    end do
end if

return
end subroutine initialise_field_configuration
