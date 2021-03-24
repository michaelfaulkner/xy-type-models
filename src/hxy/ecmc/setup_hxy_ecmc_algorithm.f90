module variables
character(100) :: output_directory
integer, parameter :: max_side = 128
integer, parameter :: max_sites = max_side * max_side
double precision, parameter :: twopi = 6.28318530717959d0
double precision, parameter :: pi = 3.14159265358979d0
double precision, parameter :: epsilon = 0.00000000001
integer :: pos_x(max_sites),neg_x(max_sites),pos_y(max_sites),neg_y(max_sites),v(max_sites)
integer :: side, sites, Tsteps, therm_sweeps, measurements, max_autocorr_time, twist, accept_twist, nmax, Nevents
integer :: calculate_external_minimising_twist_field
integer :: no_of_external_twists_to_minimise_potential_x, no_of_external_twists_to_minimise_potential_y
double precision :: theta(max_sites), top_x(max_sites), top_y(max_sites)
double precision :: sum_of_squared_electric_field_x, sum_of_squared_electric_field_y
double precision :: volume, length, T, beta, Tmin, Tmax, chainlength, maxchainlength
end module variables


! **************************************
! READ IN INPUT HYPERPARAMETERS/CONSTANTS
! **************************************

subroutine input(seed, start)
use variables
implicit none
integer seed, start

read(1, *) output_directory
read(1, *) side
read(1, *) therm_sweeps
read(1, *) measurements
read(1, *) Tmin
read(1, *) Tmax
read(1, *) Tsteps
read(1, *) start
read(1, *) twist
read(1, *) nmax
read(1, *) calculate_external_minimising_twist_field
read(1, *) seed

sites = side * side
volume = float(sites)
length = float(side)
maxchainlength = volume * twopi / 2.0 ! 2PI TO TRANSFORM TO SPIN SPACE; 1/2 TO REFLECT MANON MICHEL

if (side > max_side) then
    write(6, *) 'Linear lattice length exceeds maximum: change the maximum in the common file.'
    stop
end if

return
end subroutine input
