! **************************************
! SET VARIABLES
! **************************************

module variables
  character(100) :: output_directory
  real*8, parameter :: twopi = 6.28318530718
  real*8, parameter :: epsilon = 10.0 ** (-6)
  integer max_side,max_sites,max_sweeps
  parameter (max_side = 128)
  parameter (max_sites = max_side * max_side)
  integer*8 Nevents
  integer pos_x(max_sites),neg_x(max_sites),pos_y(max_sites),neg_y(max_sites),v(max_sites)
  integer side, sites, Tsteps, therm_sweeps, measurements, max_autocorr_time, twist, accept_twist, nmax
  integer calculate_external_minimising_twist_field
  integer  no_of_external_twists_to_minimise_potential_x, no_of_external_twists_to_minimise_potential_y
  real*8  theta(max_sites), top_x(max_sites), top_y(max_sites)
  real*8  sum_of_squared_electric_field_x, sum_of_squared_electric_field_y
  real*8  volume, length, T, beta, Tmin, Tmax, chainlength, maxchainlength
end module variables


! **************************************
! READ IN INPUT HYPERPARAMETERS/CONSTANTS
! **************************************

subroutine input(seed, start)
use variables
implicit none
integer seed,start

read(1,*) side
read(1,*) therm_sweeps
read(1,*) measurements
read(1,*) Tmin
read(1,*) Tmax
read(1,*) Tsteps
read(1,*) start
read(1,*) twist
read(1,*) nmax
read(1,*) calculate_external_minimising_twist_field
read(1,*) seed
read(1, *) output_directory

sites = side * side
volume = float(sites)
length = float(side)
maxchainlength = volume * twopi / 2.0 ! 2PI TO TRANSFORM TO SPIN SPACE; 1/2 TO REFLECT MANON MICHEL

if (side > max_side) then
    write(6,*) 'Linear lattice length exceeds maximum: change the maximum in the common file.'
    stop
end if

return
end subroutine input
