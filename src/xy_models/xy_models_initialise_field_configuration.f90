subroutine initialise_field_configuration
use variables
implicit none
integer :: i

if (randomise_initial_field_configuration) then
    do i = 1, no_of_sites
        theta(i) = twopi * rand()
    end do
else
    do i = 1, no_of_sites
        theta(i) = 0.0d0
    end do
end if

return
end subroutine initialise_field_configuration
