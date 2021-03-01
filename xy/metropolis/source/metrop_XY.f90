! **************************************
! Metropolis sampling from XY distribution
! **************************************

! **************************************
! Main program
! **************************************

program metrop_XY

  use variables
  implicit none
  integer i,j,seed,start
  real*8 Tincr

  open (unit = 1, file = 'initial.in')                                                   ! CALL INPUT HYPERPARAMTERS AND SET UP LATTICE, SEED, ETC.
  call input(seed,start)
  call PBC
  call randinit(seed)
  write(6,*) rand(seed)

  T = Tmin
  if (Tsteps.eq.0) then
     Tincr = 0.0
  else
     Tincr = (Tmax - Tmin) / Tsteps
  end if

  do i = 0, Tsteps
     write(6,*) T
     beta = 1 / T
     if (T.eq.Tmin) then
        call initial_spins(start)                                                   ! INITIALIZE SYSTEM
     end if

     do j = 0,therm_sweeps - 1                                                      ! THERMALIZE SYSTEM AT NEW TEMPERATURE
        call markov_chain_XY
     end do

     call initial_measure                                                           ! SET ALL MEASUREMENT DATA TO ZERO

     do j = 0, measurements - 1                                                     ! RUN EVENT CHAIN UNTIL MAXIMUM EVENT-CHAIN LENGTH
        call markov_chain_XY                                                        ! IS REACHED A measurements NUMBER OF TIMES
        if (twist .eq. 1) then                                                      ! DIRECT TSF PROPOSAL VIA GLOBAL TWIST, WHICH ENSURES ERGODICITY (METROPOLIS MCMC)
           call global_twist_XY
        end if
        call measure                                                                ! MEASURE THE STATE OF THE SYSTEM
     end do

     call output_acceptance_rates

     T = T + Tincr                                                                  ! SET NEW TEMPERATURE
     proposalInterval = proposalInterval + deltaProposalInterval
  end do

end program metrop_XY  

! **************************************
! METROPOLIS MARKOV-CHAIN SUBROUTINE
! **************************************

subroutine markov_chain_XY
  use variables
  implicit none
  integer n, i
  real*8 thetaOld, thetaNew, deltaTheta, Uold, Unew, deltaU
  real*8 thetaPos_x, thetaNeg_x, thetaPos_y, thetaNeg_y

   do n = 0, sites - 1
      i = int(rand() * sites)

      thetaOld = theta(i)
      deltaTheta = 2. * proposalInterval * (rand() - 0.5)
      thetaNew = thetaOld + deltaTheta

      ! CALL OTHER RELEVANT SPINS

      thetaPos_x = theta(pos_x(i))
      thetaNeg_x = theta(neg_x(i))
      thetaPos_y = theta(pos_y(i))
      thetaNeg_y = theta(neg_y(i))

      ! METROPOLIS FILTER

      Uold = - cos(thetaPos_x - thetaOld) - cos(thetaPos_y - thetaOld) - cos(thetaOld - thetaNeg_x) - cos(thetaOld - thetaNeg_y)
      Unew = - cos(thetaPos_x - thetaNew) - cos(thetaPos_y - thetaNew) - cos(thetaNew - thetaNeg_x) - cos(thetaNew - thetaNeg_y)
      deltaU = Unew - Uold

      if ((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) then
         if (thetaNew .le. -pi) then
            thetaNew = thetaNew + twopi
         else if (thetaNew .gt. pi) then
            thetaNew = thetaNew - twopi
         end if
         theta(i) = thetaNew
         accept = accept + 1
      end if
   end do
  
  return
end subroutine markov_chain_XY
