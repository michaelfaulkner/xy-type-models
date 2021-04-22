module variables
character(100) :: output_directory, algorithm_name
logical :: use_external_global_moves, randomise_initial_field_configuration, calculate_external_minimising_twist_field
integer, allocatable, dimension(:) :: pos_x, neg_x, pos_y, neg_y
integer :: integer_lattice_length, no_of_sites, no_of_temperature_increments, no_of_equilibration_sweeps
integer :: no_of_observations, no_of_events, no_of_accepted_external_global_moves, vacuum_permittivity_sum_cutoff
integer :: no_of_external_twists_to_minimise_potential_x, no_of_external_twists_to_minimise_potential_y
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, parameter :: pi_squared = 9.86960440108936d0
double precision, parameter :: pi_squared_over_two = 4.93480220054468d0
double precision, allocatable, dimension(:) :: theta, emergent_field_x, emergent_field_y
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: spin_space_distance_between_observations
double precision :: sum_of_squared_electric_field_x, sum_of_squared_electric_field_y
end module variables


subroutine read_in_config_file
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
allocate(theta(no_of_sites), emergent_field_x(no_of_sites), emergent_field_y(no_of_sites))
allocate(pos_x(no_of_sites), pos_y(no_of_sites), neg_x(no_of_sites), neg_y(no_of_sites))
spin_space_distance_between_observations = dfloat(no_of_sites) * pi

return
end subroutine read_in_config_file
