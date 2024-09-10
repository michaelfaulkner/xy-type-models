subroutine process_sample(sample_index)
use variables
implicit none
integer :: i, sample_index
double precision :: potential, non_normalised_magnetisation(2), magnetic_norm_squared

if (measure_magnetisation) then
    non_normalised_magnetisation = (/ 0.0d0, 0.0d0 /)
    do i = 1, no_of_sites
        non_normalised_magnetisation(1) = non_normalised_magnetisation(1) + cos(spin_field(i))
        non_normalised_magnetisation(2) = non_normalised_magnetisation(2) + sin(spin_field(i))
    end do
    if (print_samples) then
        write(30, 100) non_normalised_magnetisation(1), non_normalised_magnetisation(2)
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        ! nb, we divide by no_of_sites at this point to avoid infinities at large system size and long timescale
        magnetic_norm_squared = non_normalised_magnetisation(1) ** 2 / dfloat(no_of_sites ** 2) + &
                                non_normalised_magnetisation(2) ** 2 / dfloat(no_of_sites ** 2)
        magnetic_norm_sum = magnetic_norm_sum + magnetic_norm_squared ** 0.5
        magnetic_norm_squared_sum = magnetic_norm_squared_sum + magnetic_norm_squared
        magnetic_norm_quartic_sum = magnetic_norm_quartic_sum + magnetic_norm_squared ** 2
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
    if (sample_index >= no_of_equilibration_sweeps) then
        potential_sum = potential_sum + potential
        potential_squared_sum = potential_squared_sum + potential ** 2
    end if
end if

100 format(ES25.14, ",", ES25.14)
200 format(ES25.14)

return
end subroutine process_sample
