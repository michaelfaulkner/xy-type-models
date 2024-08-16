subroutine read_config_file
use variables
implicit none

read(10, *) algorithm_name
read(10, *) output_directory
read(10, *) integer_lattice_length
read(10, *) no_of_equilibration_sweeps
read(10, *) no_of_samples
read(10, *) initial_temperature
read(10, *) final_temperature
read(10, *) no_of_temperature_increments
read(10, *) vacuum_permittivity_sum_cutoff
read(10, *) randomise_initial_field_configuration
read(10, *) use_external_global_moves
read(10, *) measure_magnetisation
read(10, *) measure_helicity
read(10, *) measure_potential
read(10, *) measure_hot_twist_relaxations
read(10, *) measure_external_global_moves
read(10, *) print_samples

if (algorithm_name /= 'hxy-ecmc') then
   write(6, *) 'ConfigurationError: the value of algorithm_name does not equal hxy-ecmc.'
   stop
end if

measure_twist_relaxations = .false.
if (.not.(use_external_global_moves)) then
    measure_external_global_moves = .false.
end if
if (measure_helicity) then
    measure_emergent_field = .true.
else
    measure_emergent_field = .false.
end if
no_of_sites = integer_lattice_length * integer_lattice_length
allocate(spin_field(no_of_sites), emergent_field(no_of_sites, 2))
allocate(get_north_neighbour(no_of_sites), get_south_neighbour(no_of_sites))
allocate(get_east_neighbour(no_of_sites), get_west_neighbour(no_of_sites))
allocate(get_up_neighbour(no_of_sites), get_down_neighbour(no_of_sites))

return
end subroutine read_config_file
