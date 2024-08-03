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
read(10, *) charge_hop_proportion
read(10, *) use_external_global_moves
read(10, *) measure_electric_field_sum
read(10, *) measure_potential
read(10, *) measure_external_global_moves
read(10, *) print_samples

if ((algorithm_name /= 'elementary-electrolyte').and.(algorithm_name /= 'multivalued-electrolyte')) then
   write(6, *) 'ConfigurationError: the value of algorithm_name does not equal either elementary-electrolyte or &
                    multivalued-electrolyte.'
   stop
end if

if (.not.(use_external_global_moves)) then
    measure_external_global_moves = .false.
end if
no_of_sites = integer_lattice_length * integer_lattice_length
allocate(electric_field(no_of_sites, 2))
allocate(get_north_neighbour(no_of_sites), get_south_neighbour(no_of_sites))
allocate(get_east_neighbour(no_of_sites), get_west_neighbour(no_of_sites))
allocate(get_up_neighbour(no_of_sites), get_down_neighbour(no_of_sites), array_of_sites(no_of_sites))
allocate(charge_configuration(no_of_sites))
do i = 1, no_of_sites
    array_of_sites(i) = i
end do
charge_hop_proportion_over_two = 0.5d0 * charge_hop_proportion

return
end subroutine read_config_file
