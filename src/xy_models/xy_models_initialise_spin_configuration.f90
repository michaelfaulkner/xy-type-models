subroutine initialise_field_configuration(start)
use variables
implicit none
integer :: i, start

if (start == 0) then
    do i = 1, sites
        theta(i) = 0.0d0
    end do
else
    do i = 1, sites
        theta(i) = twopi * (rand() - 0.5d0)
    end do
end if

return
end subroutine initialise_field_configuration
