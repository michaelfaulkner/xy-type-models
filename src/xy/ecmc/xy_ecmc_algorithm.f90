program xy_ecmc_algorithm
use variables
implicit none
character(100) :: config_file
integer :: i, j, seed, start
double precision :: magnitude_of_temperature_increments

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
write(6, '(A, F16.14)') 'Initial random number = ', rand(seed)

T = Tmin
if (Tsteps == 0) then
    magnitude_of_temperature_increments = 0.0
else
    magnitude_of_temperature_increments = (Tmax - Tmin) / Tsteps
end if

do i = 0, Tsteps
    write(6, '(A, ES8.2)') 'Temperature = ', T
    beta = 1.0 / T

    do j = 1, therm_sweeps
        call event_chain_XY
    end do

    Nevents = 0
    accept_twist = 0
    call create_sample_files

    do j = 1, measurements
        call event_chain_XY
        if (twist == 1) then
            call global_twist_XY
        end if
        call measure
    end do

    call output_Nevents
    T = T + magnitude_of_temperature_increments
end do
end program xy_ecmc_algorithm


! **************************************
! EVENT-CHAIN SUBROUTINE
! **************************************

subroutine event_chain_XY
  use variables
  implicit none
  integer :: i,lift,veto
  integer, dimension (1:4) :: neighbours
  double precision :: Estar,distanceToNextEvent,deltaEinitial,deltaEexit,distanceToGo
  double precision :: thetalift,thetafix,deltaThetaInitial,Ntours,deltaThetaExit,distance

  lift = int(volume * rand())                                                       ! PICK A SITE AT RANDOM FOR INITIAL LIFTING VARIABLE (FROM {0, 1, ..., N-1})
  distanceToGo = maxchainlength                                                     ! SET REMAINING EVENT-CHAIN LENGTH AS THE MAX CHAIN LENGTH

  do                                                                                ! ITERATE UNTIL THE EVENT-CHAIN LENGTH EXCEEDS THE MAXIMUM PERMITTED: THEN WE EXIT SUBROUTINE, MEASURE AND RESAMPLE THE LIFTING SPIN/PARAMETER
     thetalift = theta(lift)

     neighbours(1) = neg_x(lift)
     neighbours(2) = pos_x(lift)
     neighbours(3) = neg_y(lift)
     neighbours(4) = pos_y(lift)

     distanceToNextEvent = 10d10                                                    ! RESET DISTANCE COUNTER FOR NEW SUB-CHAIN (WHERE SUB-CHAINS CONNECT EVENTS WITHIN AN EVENT CHAIN)

     do i = 1, 4                                                                       ! ITERATE OVER lift'S NEIGHBOURING SPINS/PARAMETERS

        thetafix = theta(neighbours(i))                                             ! THE VALUE OF EACH NEIGHBOURING SPIN/PARAMETER
        deltaThetaInitial = modulo(thetalift - thetafix + pi, twopi) - pi           ! CALCULATE SPIN/PARAMETER DIFFERENCE MODULO WRT (-PI,PI] (FORTRAN MOD FUNCTION FAILS FOR NEGATIVE VALUES HENCE THE IF LOOP)

        Estar = 1.0 - rand()                                                        ! DRAW UNIFORM CONT. RANDOM VARIABLE FROM (0,1] (TRANSFORM TO AVOIDS DIVERGENCES AT -ln(0))
        Estar = - T * log(Estar)

        if (deltaThetaInitial > 0) then                                            ! IF DERIVATIVE OF CONDITIONAL INT. POTENTIAL IS +VE
           deltaEinitial = 1 - cos(deltaThetaInitial)                               ! INITIAL CONDITIONAL INT. POTENTIAL
           Ntours = int((deltaEinitial + Estar) / 2.0d0)                              ! NO. COMPLETE ROTATIONS INDUCED THROUGH SPIN SPACE
           deltaEexit = (Ntours + 1) * 2.0d0 - (deltaEinitial + Estar)                ! CONDITIONAL-INT.-POTENTIAL DIFF. BETWEEN TOP OF NEXT MAX. AND MOD(deltaEinitial + Estar,2.0)
           deltaThetaExit = acos(1.0d0 - deltaEexit)                                  ! EQUIV. TO ABOVE LINE IN SPIN/PARAMETER SPACE
           distance = (Ntours + 0.5d0) * twopi - (deltaThetaInitial + deltaThetaExit) ! TOTAL LENGTH OF PATH COVERED IN PARAMETER SPACE
        else                                                                        ! IF DERIVATIVE OF CONDITIONAL INT. POTENTIAL IS -VE, GO TO THE BOTTOM OF THE WELL
           Ntours = int(Estar / 2.0d0)
           deltaEexit = (Ntours + 1) * 2.0d0 - Estar
           deltaThetaExit = acos(1.0d0 - deltaEexit)
           distance = (Ntours + 0.5d0) * twopi - (deltaThetaInitial + deltaThetaExit)
        end if

        if (distanceToNextEvent > distance) then                                   ! IF LENGTH COVERED IS SHORTEST SO FAR
           distanceToNextEvent = distance
           veto = neighbours(i)                                                     ! UPDATE veto SPIN/PARAMETER
        end if

     end do

     if (distanceToNextEvent > distanceToGo) then                                  ! IF MAX. EVENT-CHAIN LENGTH HAS BEEN EXCEEDED
        theta(lift) = mod(thetalift + distanceToGo, twopi)            ! JUST ADD THE REMAINING LENGTH FROM TOTAL ALLOWED CHAIN LENGTH
        exit                                                                        ! EXIT SUBROUTINE AS MAX. EVENT-CHAIN LENGTH HAS BEEN EXCEEDED: NOW MEASURE THE SYSTEM AND RESAMPLE LIFTING SPIN/PARAMETER
     else
        theta(lift) = mod(thetalift + distanceToNextEvent, twopi)      ! FINAL VALUE OF LIFTING VARIABLE IF MAX. LENGTH HASN'T BEEN EXCEEDED
        distanceToGo = distanceToGo - distanceToNextEvent
        lift = veto                                                                 ! UPDATE THE LIFTING SPIN/PARAMETER TO THAT WHICH VETOED THE CURRENT MOVE
        Nevents = Nevents + 1
     end if
  end do

  return
end subroutine event_chain_XY

! **************************************
! OUTPUT Nevents
! **************************************

subroutine output_Nevents
use variables
implicit none
character(100), parameter :: temperature_string="/temp_eq_"
character(100) :: filename

write (filename, '(A, F4.2, "//number_of_events.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit = 300, file = filename)
write(300, 100) Nevents
if (twist == 1) then
    write(300, 100) accept_twist
end if
close(300)

100 format(I20)

return
end subroutine output_Nevents
