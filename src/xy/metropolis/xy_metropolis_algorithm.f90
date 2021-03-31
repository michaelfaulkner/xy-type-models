program xy_metropolis_algorithm
use variables
implicit none
character(100) :: config_file
integer i, j, seed, start
double precision Tincr

! verify that the something has been parsed to the exectuable
if (command_argument_count() /= 1) then
    write(6, *) 'Error: parse configuration file to executable'
    stop
end if
! read in config file
call get_command_argument(1, config_file)
open (unit=1, file=config_file)

call input(seed, start)
call setup_periodic_boundaries
call initialise_spin_configuration(start)
call randinit(seed)
write(6, *) rand(seed)

T = Tmin
if (Tsteps == 0) then
    Tincr = 0.0
else
    Tincr = (Tmax - Tmin) / Tsteps
end if

do i = 0, Tsteps
    write(6, *) T
    beta = 1.0 / T

    do j = 1, therm_sweeps
        call markov_chain_XY
    end do

    accept = 0
    accept_twist = 0
    call create_sample_files

    do j = 1, measurements
        call markov_chain_XY
        if (twist == 1) then
            call global_twist_XY
        end if
        call measure
    end do

    call output_acceptance_rates
    T = T + Tincr
    proposalInterval = proposalInterval + deltaProposalInterval
end do
end program xy_metropolis_algorithm

! **************************************
! METROPOLIS MARKOV-CHAIN SUBROUTINE
! **************************************

subroutine markov_chain_XY
  use variables
  implicit none
  integer :: n, i
  double precision :: thetaOld, thetaNew, deltaTheta, Uold, Unew, deltaU
  double precision :: thetaPos_x, thetaNeg_x, thetaPos_y, thetaNeg_y

   do n = 1, sites
      i = int(rand() * sites)

      thetaOld = theta(i)
      deltaTheta = 2.0d0 * proposalInterval * (rand() - 0.5d0)
      thetaNew = mod(thetaOld + deltaTheta, twopi)

      ! CALL OTHER RELEVANT SPINS

      thetaPos_x = theta(pos_x(i))
      thetaNeg_x = theta(neg_x(i))
      thetaPos_y = theta(pos_y(i))
      thetaNeg_y = theta(neg_y(i))

      ! METROPOLIS FILTER

      Uold = - cos(thetaPos_x - thetaOld) - cos(thetaPos_y - thetaOld) - cos(thetaOld - thetaNeg_x) - cos(thetaOld - thetaNeg_y)
      Unew = - cos(thetaPos_x - thetaNew) - cos(thetaPos_y - thetaNew) - cos(thetaNew - thetaNeg_x) - cos(thetaNew - thetaNeg_y)
      deltaU = Unew - Uold

      if ((deltaU < 0.0d0) .or. (rand() < exp(- beta * deltaU))) then
         theta(i) = thetaNew
         accept = accept + 1
      end if
   end do
  
  return
end subroutine markov_chain_XY


subroutine output_acceptance_rates
use variables
implicit none
character(100), parameter :: temperature_string="/temp_eq_"
character(100) :: filename
double precision :: acceptanceRate, twistAcceptanceRate

acceptanceRate = float(accept) / (measurements * volume)
twistAcceptanceRate = float(accept_twist) / (measurements * volume)

write (filename, '(A, F4.2, "//acceptance_rates.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit = 300, file = filename)

write(300, 100) acceptanceRate
if (twist .eq. 1) then
    write(300, 100) twistAcceptanceRate
end if
close(300)

100 format(F16.8)

return
end subroutine output_acceptance_rates
