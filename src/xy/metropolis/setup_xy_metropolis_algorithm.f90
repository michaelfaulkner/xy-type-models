! **************************************
! SET VARIABLES
! **************************************

module variables
  real*8, parameter :: twopi = 6.28318530718
  real*8, parameter :: pi = 3.14159265359
  real*8, parameter :: epsilon = 10.0 ** (-6)
  integer max_side, max_sites
  parameter (max_side = 128)
  parameter (max_sites = max_side * max_side)
  integer pos_x(max_sites), neg_x(max_sites), pos_y(max_sites), neg_y(max_sites), v(max_sites)
  integer side, sites, Tsteps, therm_sweeps, measurements, twist, accept, accept_twist
  real*8  theta(max_sites), top_x(max_sites), top_y(max_sites)
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
  read(1,*) therm_sweeps
  read(1,*) measurements
  read(1,*) Tmin
  read(1,*) Tmax
  read(1,*) Tsteps
  read(1,*) proposalInterval
  read(1,*) deltaProposalInterval
  read(1,*) start
  read(1,*) twist
  read(1,*) seed
  
  sites = side * side
  volume = float(sites)

  if (side.gt.max_side) then
     write(6,*) 'Linear lattice length exceeds maximum: change the maximum in module variables.'
     stop
  end if

  return
end subroutine input
