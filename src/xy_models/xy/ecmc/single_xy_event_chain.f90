subroutine single_event_chain
use variables
implicit none
integer :: i, active_spin_index, vetoeing_spin_index
integer, dimension (1:4) :: neighbouring_spin_indices
double precision :: Estar, distanceToNextEvent, deltaEinitial, deltaEexit, distanceToGo, active_spin, thetafix
double precision :: deltaThetaInitial, deltaThetaExit, distance, Ntours

active_spin_index = int(volume * rand())                                                       ! PICK A SITE AT RANDOM FOR INITIAL LIFTING VARIABLE (FROM {0, 1, ..., N-1})
distanceToGo = maxchainlength                                                     ! SET REMAINING EVENT-CHAIN LENGTH AS THE MAX CHAIN LENGTH

do                                                                                ! ITERATE UNTIL THE EVENT-CHAIN LENGTH EXCEEDS THE MAXIMUM PERMITTED: THEN WE EXIT SUBROUTINE, MEASURE AND RESAMPLE THE LIFTING SPIN/PARAMETER
    active_spin = theta(active_spin_index)
    neighbouring_spin_indices(1) = neg_x(active_spin_index)
    neighbouring_spin_indices(2) = pos_x(active_spin_index)
    neighbouring_spin_indices(3) = neg_y(active_spin_index)
    neighbouring_spin_indices(4) = pos_y(active_spin_index)

    distanceToNextEvent = 10d10                                             ! RESET DISTANCE COUNTER FOR NEW SUB-CHAIN (WHERE SUB-CHAINS CONNECT EVENTS WITHIN AN EVENT CHAIN)

    do i = 1, 4                                                                       ! ITERATE OVER lift'S NEIGHBOURING SPINS/PARAMETERS
        thetafix = theta(neighbouring_spin_indices(i))                                             ! THE VALUE OF EACH NEIGHBOURING SPIN/PARAMETER
        deltaThetaInitial = modulo(active_spin - thetafix + pi, twopi) - pi           ! CALCULATE SPIN/PARAMETER DIFFERENCE MODULO WRT (-PI,PI] (FORTRAN MOD FUNCTION FAILS FOR NEGATIVE VALUES HENCE THE IF LOOP)

        Estar = 1.0 - rand()                                                        ! DRAW UNIFORM CONT. RANDOM VARIABLE FROM (0,1] (TRANSFORM TO AVOIDS DIVERGENCES AT -ln(0))
        Estar = - temperature * log(Estar)

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
           vetoeing_spin_index = neighbouring_spin_indices(i)                                                     ! UPDATE veto SPIN/PARAMETER
        end if

    end do

    if (distanceToNextEvent > distanceToGo) then                                  ! IF MAX. EVENT-CHAIN LENGTH HAS BEEN EXCEEDED
        theta(active_spin_index) = mod(active_spin + distanceToGo, twopi)            ! JUST ADD THE REMAINING LENGTH FROM TOTAL ALLOWED CHAIN LENGTH
        exit                                                                        ! EXIT SUBROUTINE AS MAX. EVENT-CHAIN LENGTH HAS BEEN EXCEEDED: NOW MEASURE THE SYSTEM AND RESAMPLE LIFTING SPIN/PARAMETER
    else
        theta(active_spin_index) = mod(active_spin + distanceToNextEvent, twopi)      ! FINAL VALUE OF LIFTING VARIABLE IF MAX. LENGTH HASN'T BEEN EXCEEDED
        distanceToGo = distanceToGo - distanceToNextEvent
        active_spin_index = vetoeing_spin_index                                                                 ! UPDATE THE LIFTING SPIN/PARAMETER TO THAT WHICH VETOED THE CURRENT MOVE
        no_of_events = no_of_events + 1
    end if
end do

return
end subroutine single_event_chain
