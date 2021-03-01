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
! INITIAL SPIN CONFIGURATION
! **************************************

subroutine initial_spins(start)
  use variables
  implicit none
  integer start,i,j
  if (start.eq.0) then
     do i=0,sites-1
        theta(i) = 0.0
     end do
  else
     do i=0,sites-1
        theta(i) = twopi * rand()
     end do
  end if
  return
end subroutine initial_spins

! **************************************
! SET CHAIN LENGTH AND Nevents TO 0; CREATE NEW DIRECTORY AND FILES FOR NEW TEMP
! **************************************

subroutine initial_measure
  use variables
  implicit none
  character(100) filename
  character(100) directory

  accept = 0
  accept_twist = 0                                                                ! SET TOTAL NUMBER OF GLOBAL TWISTS TO ZERO

  write (directory, '( "temp_eq_", F4.2)' ) T                                     ! OPENS NEW DIRECTORY IN WHICH TO SAVE THE MARKOV CHAIN FOR THE CURRENT TEMPERATURE
  call system ( 'mkdir -p ' // directory )

  write (filename, '( "temp_eq_", F4.2,"//magn_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=10,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//magn_x_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=11,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//magn_y_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=12,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//cos_top_x_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=13,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//cos_top_y_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=14,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//sin_top_x_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=15,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//sin_top_y_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=16,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//Ebar_x_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=17,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//Ebar_y_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=18,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//vort_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=19,file=filename)
  write (filename, '( "temp_eq_", F4.2,"//potential_samples_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit = 20, file = filename)

  return
end subroutine initial_measure
