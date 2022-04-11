subroutine get_and_print_observation
use variables
implicit none
integer :: i
double precision :: potential, non_normalised_magnetisation(2)

if (measure_magnetisation) then
    non_normalised_magnetisation = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        non_normalised_magnetisation(1) = non_normalised_magnetisation(1) + cos(spin_field(i))
        non_normalised_magnetisation(2) = non_normalised_magnetisation(2) + sin(spin_field(i))
    end do
    write(30, 100) non_normalised_magnetisation(1), non_normalised_magnetisation(2)
end if

if (measure_potential) then
    potential = 0.0d0
    do i = 1, no_of_sites
        potential = potential - cos(spin_field(get_east_neighbour(i)) - spin_field(i)) - &
                                cos(spin_field(get_north_neighbour(i)) - spin_field(i)) - &
                                cos(spin_field(get_up_neighbour(i)) - spin_field(i))
    end do
    write(40, 200) potential
end if

100 format(ES25.14, ",", ES25.14)
200 format(ES25.14)

return
end subroutine get_and_print_observation
