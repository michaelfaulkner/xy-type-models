subroutine read_config_file
use variables
implicit none
integer :: i

read(10, *) algorithm_name
read(10, *) output_directory
read(10, *) integer_lattice_length
read(10, *) no_of_equilibration_sweeps
read(10, *) no_of_samples
read(10, *) initial_temperature
read(10, *) final_temperature
read(10, *) no_of_temperature_increments
read(10, *) width_of_proposal_interval
read(10, *) target_acceptance_rate_of_field_rotations
read(10, *) always_cold_start
read(10, *) always_hot_start
read(10, *) use_external_global_moves
read(10, *) measure_magnetisation
read(10, *) measure_potential
read(10, *) print_samples

if (algorithm_name /= '3dxy-gaussian-noise-metropolis') then
   write(6, *) 'ConfigurationError: the value of algorithm_name does not equal 3dxy-gaussian-noise-metropolis.'
   stop
end if

if ((always_cold_start).and.(always_hot_start)) then
    write(6, *) 'ConfigurationError: the value of always_cold_start and always_hot_start are both equal to true.'
   stop
end if

use_external_global_moves = .false.
no_of_sites = integer_lattice_length * integer_lattice_length * integer_lattice_length
allocate(spin_field(no_of_sites))
allocate(get_north_neighbour(no_of_sites), get_south_neighbour(no_of_sites))
allocate(get_east_neighbour(no_of_sites), get_west_neighbour(no_of_sites))
allocate(get_up_neighbour(no_of_sites), get_down_neighbour(no_of_sites), array_of_sites(no_of_sites))
do i = 1, no_of_sites
    array_of_sites(i) = i
end do
! charge_hop_proportion is a dummy variable that allows the electrolyte and h/xy models to use the same
! step-size adaptor in the same metropolis_algorithm.f90 file
charge_hop_proportion = 0.0d0

return
end subroutine read_config_file
