module variables
character(100) :: output_directory, algorithm_name
logical :: use_external_global_moves, randomise_initial_field_configuration
integer, parameter :: max_side = 128
integer, parameter :: max_sites = max_side * max_side
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
integer :: pos_x(max_sites), neg_x(max_sites), pos_y(max_sites), neg_y(max_sites), v(max_sites)
integer :: side, sites, no_of_temperature_increments, therm_sweeps, measurements, max_autocorr_time
integer :: no_of_events, no_of_accepted_external_global_moves
double precision :: theta(max_sites), volume, length, beta, temperature, initial_temperature, final_temperature
double precision :: magnitude_of_temperature_increments, spin_space_distance_between_observations
end module variables


subroutine read_in_config_file
use variables
implicit none

read(10, *) algorithm_name
read(10, *) output_directory
read(10, *) side
read(10, *) therm_sweeps
read(10, *) measurements
read(10, *) initial_temperature
read(10, *) final_temperature
read(10, *) no_of_temperature_increments
read(10, *) randomise_initial_field_configuration
read(10, *) use_external_global_moves

if (algorithm_name /= 'xy-ecmc') then
   write(6, *) 'ConfigurationFileError: the value of algorithm_name does not equal xy-ecmc.'
   stop
end if

sites = side * side
length = float(side)
volume = float(sites)
spin_space_distance_between_observations = volume * pi

if (side > max_side) then
   write(6, *) 'Linear lattice length exceeds maximum: change the maximum in module variables.'
   stop
end if

return
end subroutine read_in_config_file
