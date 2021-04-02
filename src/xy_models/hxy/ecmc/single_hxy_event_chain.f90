subroutine single_event_chain
  use variables
  implicit none
  integer :: i, lift, veto
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

     do i = 1, 4                                                                    ! ITERATE OVER lift'S NEIGHBOURING SPINS/PARAMETERS
        thetafix = theta(neighbours(i))                                             ! THE VALUE OF EACH NEIGHBOURING SPIN/PARAMETER
        deltaThetaInitial = modulo(thetalift - thetafix + pi, twopi) - pi           ! CALCULATE SPIN/PARAMETER DIFFERENCE MODULO WRT (-PI,PI] (FORTRAN MOD FUNCTION FAILS FOR NEGATIVE VALUES HENCE THE IF LOOP)
        Estar = 1.0d0 - rand()                                                        ! DRAW UNIFORM CONT. RANDOM VARIABLE FROM (0,1] (TRANSFORM TO AVOIDS DIVERGENCES AT -ln(0))
        Estar = - temperature * log(Estar)

        if (deltaThetaInitial > 0.0d0) then                                                     ! IF DERIVATIVE OF CONDITIONAL INT. POTENTIAL IS +VE
           deltaEinitial = 0.5d0 * deltaThetaInitial * deltaThetaInitial                        ! INITIAL CONDITIONAL INT. POTENTIAL
           Ntours = dble(int((deltaEinitial + Estar) / (twopi * twopi / 8.0d0)))                ! NO. COMPLETE ROTATIONS INDUCED THROUGH SPIN SPACE
           deltaEexit = (Ntours + 1.0d0) * twopi * twopi / 8.0d0 - (deltaEinitial + Estar)      ! CONDITIONAL-INT.-POTENTIAL DIFF. BETWEEN TOP OF NEXT MAX. AND MOD(deltaEinitial + Estar,2.0)
           deltaThetaExit = pi - sqrt(twopi * twopi / 4.0d0 - 2.0d0 * deltaEexit)               ! EQUIV. TO ABOVE LINE IN SPIN/PARAMETER SPACE
           distance = (Ntours + 0.5d0) * twopi - (deltaThetaInitial + deltaThetaExit)           ! TOTAL LENGTH OF PATH COVERED IN PARAMETER SPACE
        else                                                                                    ! IF DERIVATIVE OF CONDITIONAL INT. POTENTIAL IS -VE, GO TO THE BOTTOM OF THE WELL
           Ntours = dble(int(Estar / (twopi * twopi / 8.0d0)))
           deltaEexit = (Ntours + 1.0d0) * twopi * twopi / 8.0d0 - Estar
           deltaThetaExit = twopi / 2.0d0 - sqrt(twopi * twopi / 4.0d0 - 2.0d0 * deltaEexit)
           distance = (Ntours + 0.5d0) * twopi - (deltaThetaInitial + deltaThetaExit)
        end if

        if (distanceToNextEvent > distance) then                                 ! IF LENGTH COVERED IS SHORTEST SO FAR
           distanceToNextEvent = distance
           veto = neighbours(i)                                                     ! UPDATE veto SPIN/PARAMETER
        end if

     end do

     if (distanceToNextEvent > distanceToGo) then                                  ! IF MAX. EVENT-CHAIN LENGTH HAS BEEN EXCEEDED
        theta(lift) = mod(thetalift + distanceToGo, twopi)              ! JUST ADD THE REMAINING LENGTH FROM TOTAL ALLOWED CHAIN LENGTH
        chainlength = chainlength + distanceToGo
        exit                                                                        ! EXIT SUBROUTINE AS MAX. EVENT-CHAIN LENGTH HAS BEEN EXCEEDED: NOW MEASURE THE SYSTEM AND RESAMPLE LIFTING SPIN/PARAMETER
     else
        theta(lift) = mod(thetalift + distanceToNextEvent, twopi)               ! FINAL VALUE OF LIFTING VARIABLE IF MAX. LENGTH HASN'T BEEN EXCEEDED
        chainlength = chainlength + distanceToNextEvent
        distanceToGo = distanceToGo - distanceToNextEvent
        lift = veto                                                                 ! UPDATE THE LIFTING SPIN/PARAMETER TO THAT WHICH VETOED THE CURRENT MOVE
        no_of_events = no_of_events + 1
     end if
  end do
  
  return
end subroutine single_event_chain
