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
integer pos_x(max_sites), neg_x(max_sites), pos_y(max_sites), neg_y(max_sites), rho(max_sites)
integer side, sites, Tsteps, thermSweeps, measurements, globalTSFon
integer accept_charge, accept_aux_field, accept_TSF, ratio_charge_updates, ratio_TSF_updates
real*8  Efield_x(max_sites), Efield_y(max_sites), Esum_x, Esum_y, elementaryCharge
real*8  volume, length, T, beta, Tmin, Tmax, proposalInterval, deltaProposalInterval
end module variables


! **************************************
! READ IN INPUT HYPERPARAMETERS/CONSTANTS
! **************************************

subroutine input(seed, start)
use variables
implicit none
integer seed, start

read(1,*) side
read(1,*) thermSweeps
read(1,*) measurements
read(1,*) Tmin
read(1,*) Tmax
read(1,*) Tsteps
read(1,*) proposalInterval
read(1,*) deltaProposalInterval
read(1,*) ratio_charge_updates
read(1,*) ratio_TSF_updates
read(1,*) globalTSFon
read(1,*) elementaryCharge
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


! **************************************
! SETS NEIGHBOURS WITH PERIODIC BOUNDARY CONDITIONS
! **************************************

subroutine PBC
  use variables
  implicit none
  integer i
  ! mod(i + side + 1,side) RATHER THAN mod(i + 1,side), ETC. BELOW AS mod(x,side)
  ! DOESN'T RETURN VALUES IN THE INTERVAL [0,side) FOR NEGATIVE x
  do i = 0,sites - 1
     pos_x(i) = i + mod(i + side + 1,side) - mod(i + side,side)
     neg_x(i) = i + mod(i + side - 1,side) - mod(i + side,side)
     pos_y(i) = i + (mod(int(i / side) + side + 1,side) - mod(int(i / side) + side,side)) * side
     neg_y(i) = i + (mod(int(i / side) + side - 1,side) - mod(int(i / side) + side,side)) * side
  end do
  return
end subroutine PBC

! **************************************
! INITIAL FIELD CONFIGURATION
! **************************************

subroutine initial_Efield
  use variables
  implicit none
  integer i

  do i = 0, sites - 1
     Efield_x(i) = 0.0
     Efield_y(i) = 0.0
  end do

  return
end subroutine initial_Efield

! **************************************
! CREATE NEW DIRECTORY AND FILES FOR NEW TEMP
! **************************************

subroutine initial_measure
use variables
implicit none
character(100) filename
character(100) temperature_directory
character(100), parameter :: temperature_string="/temp_eq_"

accept_charge = 0
accept_aux_field = 0
accept_TSF = 0

! OPENS NEW DIRECTORY IN WHICH TO SAVE THE MARKOV CHAIN FOR THE CURRENT TEMPERATURE
write (temperature_directory, '(A, F4.2)') trim(output_directory)//trim(temperature_string), T
call system ( 'mkdir -p ' // temperature_directory )

write (filename, '(A, F4.2, "//Esum_x_sample.dat")') trim(output_directory)//trim(temperature_string), T
open(unit=10, file=filename)
write (filename, '(A, F4.2, "//Esum_y_sample.dat")') trim(output_directory)//trim(temperature_string), T
open(unit=11, file=filename)
write (filename, '(A, F4.2, "//potential_sample.dat")') trim(output_directory)//trim(temperature_string), T
open(unit=12, file=filename)

return
end subroutine initial_measure
