module variables
character(100) :: output_directory, algorithm_name
integer, parameter :: max_side = 128
integer, parameter :: max_sites = max_side * max_side
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
integer :: pos_x(max_sites), neg_x(max_sites), pos_y(max_sites), neg_y(max_sites), v(max_sites)
integer :: start, side, sites, no_of_temperature_increments, therm_sweeps, measurements, max_autocorr_time, twist
integer :: no_of_events, no_of_accepted_external_global_moves
double precision :: theta(max_sites), top_x(max_sites), top_y(max_sites)
double precision :: volume, length, beta, temperature, initial_temperature, final_temperature
double precision :: spin_space_distance_between_observations
end module variables

! **************************************
! READ IN INPUT HYPERPARAMETERS/CONSTANTS
! **************************************

subroutine input(seed)
use variables
implicit none
integer :: seed

read(1, *) algorithm_name
read(1, *) output_directory
read(1, *) side
read(1, *) therm_sweeps
read(1, *) measurements
read(1, *) initial_temperature
read(1, *) final_temperature
read(1, *) no_of_temperature_increments
read(1, *) start
read(1, *) twist
read(1, *) seed

if (algorithm_name /= 'xy-ecmc') then
   write(6, *) 'ConfigurationFileError: the value of algorithm_name does not equal xy-ecmc.'
   stop
end if
temperature = initial_temperature
sites = side * side
volume = float(sites)
spin_space_distance_between_observations = volume * pi

if (side > max_side) then
   write(6, *) 'Linear lattice length exceeds maximum: change the maximum in module variables.'
   stop
end if

return
end subroutine input
