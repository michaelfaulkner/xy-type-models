module variables
character(100) :: output_directory
integer, parameter :: max_side = 128
integer, parameter :: max_sites = max_side * max_side
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
integer :: pos_x(max_sites), neg_x(max_sites), pos_y(max_sites), neg_y(max_sites), array_of_sites(max_sites)
integer :: side, sites, no_of_temperature_increments, therm_sweeps, measurements, twist, no_of_accepted_local_moves
integer :: no_of_accepted_external_global_moves
double precision :: theta(max_sites), top_x(max_sites), top_y(max_sites)
double precision :: volume, length, beta, temperature, initial_temperature, final_temperature
double precision :: width_of_proposal_interval, magnitude_of_proposal_interval_increments
end module variables

! **************************************
! READ IN INPUT HYPERPARAMETERS/CONSTANTS
! **************************************

subroutine input(seed, start)
use variables
implicit none
integer :: i, seed, start

read(1, *) output_directory
read(1,*) side
read(1,*) therm_sweeps
read(1,*) measurements
read(1,*) initial_temperature
read(1,*) final_temperature
read(1,*) no_of_temperature_increments
read(1,*) width_of_proposal_interval
read(1,*) magnitude_of_proposal_interval_increments
read(1,*) start
read(1,*) twist
read(1,*) seed

temperature = initial_temperature
sites = side * side
volume = float(sites)

do i = 1, sites
    array_of_sites(i) = i
end do

if (side > max_side) then
   write(6,*) 'Linear lattice length exceeds maximum: change the maximum in module variables.'
   stop
end if

return
end subroutine input
