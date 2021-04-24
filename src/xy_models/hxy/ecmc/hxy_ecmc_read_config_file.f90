subroutine read_config_file
use variables
implicit none

read(10, *) algorithm_name
read(10, *) output_directory
read(10, *) integer_lattice_length
read(10, *) no_of_equilibration_sweeps
read(10, *) no_of_observations
read(10, *) initial_temperature
read(10, *) final_temperature
read(10, *) no_of_temperature_increments
read(10, *) vacuum_permittivity_sum_cutoff
read(10, *) randomise_initial_field_configuration
read(10, *) use_external_global_moves
read(10, *) calculate_external_minimising_twist_field

if (algorithm_name /= 'hxy-ecmc') then
   write(6, *) 'ConfigurationFileError: the value of algorithm_name does not equal hxy-ecmc.'
   stop
end if

no_of_sites = integer_lattice_length * integer_lattice_length
allocate(spin_field(no_of_sites), emergent_field_x(no_of_sites), emergent_field_y(no_of_sites))
allocate(get_north_neighbour(no_of_sites), get_south_neighbour(no_of_sites))
allocate(get_east_neighbour(no_of_sites), get_west_neighbour(no_of_sites))
spin_space_distance_between_observations = dfloat(no_of_sites) * pi

return
end subroutine read_config_file