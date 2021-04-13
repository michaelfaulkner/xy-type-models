module variables
integer, allocatable, dimension(:) :: pos_x, neg_x, pos_y, neg_y, array_of_sites
integer :: start, side, sites, no_of_temperature_increments, therm_sweeps, measurements, twist
integer :: no_of_accepted_local_moves, no_of_accepted_external_global_moves, no_of_events
double precision, parameter :: twopi = 6.28318530717959d0
double precision, allocatable, dimension(:) :: theta
double precision :: volume, length, beta, temperature, initial_temperature, final_temperature
double precision :: width_of_proposal_interval, magnitude_of_proposal_interval_increments
character(100) :: output_directory, algorithm_name
end module variables

! **************************************
! READ IN INPUT HYPERPARAMETERS/CONSTANTS
! **************************************

subroutine input(seed)
use variables
implicit none
integer :: i, seed

read(1, *) algorithm_name
read(1, *) output_directory
read(1, *) side
read(1, *) therm_sweeps
read(1, *) measurements
read(1, *) initial_temperature
read(1, *) final_temperature
read(1, *) no_of_temperature_increments
read(1, *) width_of_proposal_interval
read(1, *) magnitude_of_proposal_interval_increments
read(1, *) start
read(1, *) twist
read(1, *) seed

if (algorithm_name /= 'xy-metropolis') then
   write(6, *) 'ConfigurationFileError: the value of algorithm_name does not equal xy-metropolis.'
   stop
end if
temperature = initial_temperature
sites = side * side
volume = float(sites)

allocate(theta(sites), pos_x(sites), pos_y(sites), neg_x(sites), neg_y(sites), array_of_sites(sites))

do i = 1, sites
    array_of_sites(i) = i
end do

return
end subroutine input
