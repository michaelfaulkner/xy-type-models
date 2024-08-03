subroutine process_sample(observation_index)
use variables
implicit none
integer :: i, observation_index
double precision :: potential, non_normalised_magnetisation(2), raw_magnetic_norm_squared

if (measure_magnetisation) then
    non_normalised_magnetisation = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        non_normalised_magnetisation(1) = non_normalised_magnetisation(1) + cos(spin_field(i))
        non_normalised_magnetisation(2) = non_normalised_magnetisation(2) + sin(spin_field(i))
    end do
    if (print_samples) then
        write(30, 100) non_normalised_magnetisation(1), non_normalised_magnetisation(2)
    end if
    if (observation_index >= no_of_equilibration_sweeps) then
        ! magnetic_norm_squared = raw_magnetic_norm_squared / no_of_sites ** 2
        raw_magnetic_norm_squared = non_normalised_magnetisation(1) ** 2 + non_normalised_magnetisation(2) ** 2
        raw_magnetic_norm_sum = raw_magnetic_norm_sum + raw_magnetic_norm_squared ** 0.5
        raw_magnetic_norm_squared_sum = raw_magnetic_norm_squared_sum + raw_magnetic_norm_squared
        raw_magnetic_norm_quartic_sum = raw_magnetic_norm_quartic_sum + raw_magnetic_norm_squared ** 2
    end if
end if

if (measure_potential) then
    potential = 0.0d0
    do i = 1, no_of_sites
        potential = potential - cos(spin_field(get_east_neighbour(i)) - spin_field(i)) - &
                                cos(spin_field(get_north_neighbour(i)) - spin_field(i)) - &
                                cos(spin_field(get_up_neighbour(i)) - spin_field(i))
    end do
    if (print_samples) then
        write(40, 200) potential
    end if
    if (observation_index >= no_of_equilibration_sweeps) then
        potential_sum = potential_sum + potential
        potential_squared_sum = potential_squared_sum + potential ** 2
    end if
end if

100 format(ES25.14, ",", ES25.14)
200 format(ES25.14)

return
end subroutine process_sample
