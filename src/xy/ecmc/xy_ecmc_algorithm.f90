! **************************************
! Event-chain sampling from XY distribution
! **************************************

! **************************************
! Main program
! **************************************

program xy_ecmc_algorithm

  use variables
  implicit none
  integer i,j,seed,start
  real*8 Tincr

  open (unit=1,file='initial.in')                                                   ! CALL INPUT HYPERPARAMTERS AND SET UP LATTICE, SEED, ETC.
  call input(seed,start)
  call PBC
  call randinit(seed)
  write(6,*) rand(seed)

  chainlength = 0.0                                                               ! SET TOTAL CHAIN LENGTH TO ZERO
  Nevents = 0                                                                     ! SET TOTAL NUMBER OF EVENTS TO ZERO
  accept_twist = 0                                                                ! SET TOTAL NUMBER OF GLOBAL TWISTS TO ZERO
  T = Tmin
  if (Tsteps.eq.0) then
     Tincr = 0.0
  else
     Tincr = (Tmax - Tmin) / Tsteps
  end if

  do i = 0,Tsteps
     write(6,*) T
     beta=1/T
     if (T.eq.Tmin) then
        call initial_spins(start)                                                   ! INITIALIZE SYSTEM
     end if

     do j = 0,therm_sweeps - 1                                                      ! THERMALIZE SYSTEM AT NEW TEMPERATURE
        call event_chain_XY
     end do

     call initial_measure                                                           ! SET ALL MEASUREMENT DATA TO ZERO

     do j=0,measurements - 1                                                        ! RUN EVENT CHAIN UNTIL MAXIMUM EVENT-CHAIN LENGTH
        call event_chain_XY                                                         ! IS REACHED A measurements NUMBER OF TIMES
        if (twist.eq.1) then                                                        ! DIRECT TSF PROPOSAL VIA GLOBAL TWIST, WHICH ENSURES ERGODICITY (METROPOLIS MCMC)
           call global_twist_XY
        end if
        call measure                                                                ! MEASURE THE STATE OF THE SYSTEM
     end do

     call output_Nevents

     T = T + Tincr                                                                  ! SET NEW TEMPERATURE
  end do

end program xy_ecmc_algorithm


! **************************************
! EVENT-CHAIN SUBROUTINE
! **************************************

subroutine event_chain_XY
  use variables
  implicit none
  integer i,lift,veto
  integer, dimension (0:3) :: neighbours
  real*8 Estar,distanceToNextEvent,deltaEinitial,deltaEexit,distanceToGo
  real*8 thetalift,thetafix,deltaThetaInitial,Ntours,deltaThetaExit,distance

  lift = int(volume * rand())                                                       ! PICK A SITE AT RANDOM FOR INITIAL LIFTING VARIABLE (FROM {0, 1, ..., N-1})
  distanceToGo = maxchainlength                                                     ! SET REMAINING EVENT-CHAIN LENGTH AS THE MAX CHAIN LENGTH

  do                                                                                ! ITERATE UNTIL THE EVENT-CHAIN LENGTH EXCEEDS THE MAXIMUM PERMITTED: THEN WE EXIT SUBROUTINE, MEASURE AND RESAMPLE THE LIFTING SPIN/PARAMETER
     thetalift = theta(lift)

     neighbours(0) = neg_x(lift)
     neighbours(1) = pos_x(lift)
     neighbours(2) = neg_y(lift)
     neighbours(3) = pos_y(lift)

     distanceToNextEvent = 10d10                                                    ! RESET DISTANCE COUNTER FOR NEW SUB-CHAIN (WHERE SUB-CHAINS CONNECT EVENTS WITHIN AN EVENT CHAIN)

     do i=0,3                                                                       ! ITERATE OVER lift'S NEIGHBOURING SPINS/PARAMETERS

        thetafix = theta(neighbours(i))                                             ! THE VALUE OF EACH NEIGHBOURING SPIN/PARAMETER
        deltaThetaInitial = thetalift - thetafix                                    ! CALCULATE SPIN/PARAMETER DIFFERENCE MODULO WRT (-PI,PI] (FORTRAN MOD FUNCTION FAILS FOR NEGATIVE VALUES HENCE THE IF LOOP)
        if (deltaThetaInitial.gt.0.5 * twopi) then
           deltaThetaInitial = deltaThetaInitial - twopi
        else if (deltaThetaInitial.le.-0.5 * twopi) then
           deltaThetaInitial = deltaThetaInitial + twopi
        end if

        Estar = 1.0 - rand()                                                        ! DRAW UNIFORM CONT. RANDOM VARIABLE FROM (0,1] (TRANSFORM TO AVOIDS DIVERGENCES AT -ln(0))
        Estar = - T * log(Estar)

        if (deltaThetaInitial.gt.0) then                                            ! IF DERIVATIVE OF CONDITIONAL INT. POTENTIAL IS +VE
           deltaEinitial = 1 - cos(deltaThetaInitial)                               ! INITIAL CONDITIONAL INT. POTENTIAL
           Ntours = int((deltaEinitial + Estar) / 2.0)                              ! NO. COMPLETE ROTATIONS INDUCED THROUGH SPIN SPACE
           deltaEexit = (Ntours + 1) * 2.0 - (deltaEinitial + Estar)                ! CONDITIONAL-INT.-POTENTIAL DIFF. BETWEEN TOP OF NEXT MAX. AND MOD(deltaEinitial + Estar,2.0)
           deltaThetaExit = acos(1.0 - deltaEexit)                                  ! EQUIV. TO ABOVE LINE IN SPIN/PARAMETER SPACE
           distance = (Ntours + 0.5) * twopi - (deltaThetaInitial + deltaThetaExit) ! TOTAL LENGTH OF PATH COVERED IN PARAMETER SPACE
        else                                                                        ! IF DERIVATIVE OF CONDITIONAL INT. POTENTIAL IS -VE, GO TO THE BOTTOM OF THE WELL
           Ntours = int(Estar / 2.0)
           deltaEexit = (Ntours + 1) * 2.0 - Estar
           deltaThetaExit = acos(1.0 - deltaEexit)
           distance = (Ntours + 0.5) * twopi - (deltaThetaInitial + deltaThetaExit)
        end if

        if (distanceToNextEvent.gt.distance) then                                   ! IF LENGTH COVERED IS SHORTEST SO FAR
           distanceToNextEvent = distance
           veto = neighbours(i)                                                     ! UPDATE veto SPIN/PARAMETER
        end if

     end do

     if (distanceToNextEvent.gt.distanceToGo) then                                  ! IF MAX. EVENT-CHAIN LENGTH HAS BEEN EXCEEDED
        theta(lift) = mod(thetalift + distanceToGo,twopi)                           ! JUST ADD THE REMAINING LENGTH FROM TOTAL ALLOWED CHAIN LENGTH
        chainlength = chainlength + distanceToGo
        exit                                                                        ! EXIT SUBROUTINE AS MAX. EVENT-CHAIN LENGTH HAS BEEN EXCEEDED: NOW MEASURE THE SYSTEM AND RESAMPLE LIFTING SPIN/PARAMETER
     else
        theta(lift) = mod(thetalift + distanceToNextEvent,twopi)                    ! FINAL VALUE OF LIFTING VARIABLE IF MAX. LENGTH HASN'T BEEN EXCEEDED
        chainlength = chainlength + distanceToNextEvent
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
  character(100) filename

  write (filename, '( "temp_eq_", F4.2,"//Nevents_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=20,file=filename)
  write(20,100) Nevents
  if (twist.eq.1) then
     write(20,100) accept_twist
  end if
  close(20)

100 format(I20)

  return
end subroutine output_Nevents
