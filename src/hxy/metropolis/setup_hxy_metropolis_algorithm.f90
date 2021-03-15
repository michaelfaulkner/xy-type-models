! **************************************
! SET VARIABLES
! **************************************

module variables
  character(100) :: output_directory
  real*8, parameter :: twopi = 6.28318530718
  real*8, parameter :: pi = 3.14159265359
  real*8, parameter :: epsilon = 10.0 ** (-6)
  integer max_side, max_sites
  parameter (max_side = 128)
  parameter (max_sites = max_side * max_side)
  integer pos_x(max_sites), neg_x(max_sites), pos_y(max_sites), neg_y(max_sites), v(max_sites)
  integer side, sites, Tsteps, thermSweeps, measurements, twist, accept, accept_twist, nmax
  integer calculate_external_minimising_twist_field
  integer  no_of_external_twists_to_minimise_potential_x, no_of_external_twists_to_minimise_potential_y
  real*8  theta(max_sites), top_x(max_sites), top_y(max_sites)
  real*8  sum_of_squared_electric_field_x, sum_of_squared_electric_field_y
  real*8  volume, length, T, beta, Tmin, Tmax, proposalInterval, deltaProposalInterval
end module variables


! **************************************
! READ IN INPUT HYPERPARAMETERS/CONSTANTS
! **************************************

subroutine input(seed, start)
use variables
implicit none
integer seed,start

read(1,*) side
read(1,*) thermSweeps
read(1,*) measurements
read(1,*) Tmin
read(1,*) Tmax
read(1,*) Tsteps
read(1,*) proposalInterval
read(1,*) deltaProposalInterval
read(1,*) start
read(1,*) twist
read(1,*) nmax
read(1,*) seed
read(1, *) output_directory

sites = side * side
volume = float(sites)
length = float(side)

if (side > max_side) then
    write(6,*) 'Linear lattice length exceeds maximum: change the maximum in the common file.'
    stop
end if

return
end subroutine input
