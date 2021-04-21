module variables
integer, allocatable, dimension(:) :: pos_x, neg_x, pos_y, neg_y, array_of_sites
integer :: start, side, sites, no_of_temperature_increments, therm_sweeps, measurements, twist
integer :: no_of_accepted_field_rotations, no_of_accepted_external_global_moves, no_of_events
double precision, parameter :: twopi = 6.28318530717959d0
double precision, allocatable, dimension(:) :: theta
double precision :: beta, temperature, initial_temperature, final_temperature, magnitude_of_temperature_increments
double precision :: volume, length, width_of_proposal_interval, magnitude_of_proposal_interval_increments
character(100) :: output_directory, algorithm_name
end module variables

! **************************************
! READ IN INPUT HYPERPARAMETERS/CONSTANTS
! **************************************

subroutine read_in_config_file
use variables
implicit none
integer :: i

read(10, *) algorithm_name
read(10, *) output_directory
read(10, *) side
read(10, *) therm_sweeps
read(10, *) measurements
read(10, *) initial_temperature
read(10, *) final_temperature
read(10, *) no_of_temperature_increments
read(10, *) width_of_proposal_interval
read(10, *) magnitude_of_proposal_interval_increments
read(10, *) start
read(10, *) twist

if (algorithm_name /= 'xy-metropolis') then
   write(6, *) 'ConfigurationFileError: the value of algorithm_name does not equal xy-metropolis.'
   stop
end if

sites = side * side
length = float(side)
volume = float(sites)
allocate(theta(sites), pos_x(sites), pos_y(sites), neg_x(sites), neg_y(sites), array_of_sites(sites))
do i = 1, sites
    array_of_sites(i) = i
end do

return
end subroutine read_in_config_file
