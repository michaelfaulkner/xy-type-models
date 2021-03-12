! **************************************
! Metropolis sampling from HXY distribution
! **************************************

! **************************************
! Main program
! **************************************

program metrop_HXY

  use variables
  implicit none
  integer i, j, seed, start
  real*8 Tincr

  open (unit = 1, file = 'initial.in')
  call input(seed, start)
  call PBC
  call randinit(seed)
  write(6,*) rand(seed)
  
  T = Tmin
  if (Tsteps .eq. 0) then
     Tincr = 0.0
  else
     Tincr = (Tmax - Tmin) / Tsteps
  end if
  
  do i = 0, Tsteps
     write(6,*) T
     beta=1/T
     if (T.eq.Tmin) then
        call initial_spins(start)
     end if
     
     do j = 0, thermSweeps - 1
        call markov_chain_HXY
     end do

     call initial_measure
     
     do j = 0, measurements - 1
        call markov_chain_HXY
        if (twist .eq. 1) then
           call global_twist_HXY
        end if
        call measure
     end do

     call output_acceptance_rates

     T = T + Tincr
     proposalInterval = proposalInterval + deltaProposalInterval
  end do

end program metrop_HXY


! **************************************
! METROPOLIS MARKOV-CHAIN SUBROUTINE
! **************************************

subroutine markov_chain_HXY
  use variables
  implicit none
  integer n, i
  real*8 thetaOld, thetaNew, deltaTheta, Uold, Unew, deltaU
  real*8 top1old, top2old, top3old, top4old, top1new, top2new, top3new, top4new

   do n = 0, sites - 1
      i = int(rand() * sites)

      thetaOld = theta(i)
      deltaTheta = 2. * proposalInterval * (rand() - 0.5)
      thetaNew = thetaOld + deltaTheta
      if (thetaNew .le. -pi) then
         thetaNew = thetaNew + twopi
      else if (thetaNew .gt. pi) then
         thetaNew = thetaNew - twopi
      end if

      ! CALL OLD EMERGENT FIELD

      top1old = top_x(i)
      top2old = top_y(i)
      top3old = top_x(pos_y(i))
      top4old = top_y(pos_x(i))

      ! PROPOSED EMERGENT FIELD
      
      top1new = thetanew - theta(neg_y(i))
      if (top1new .gt. pi) then
         top1new = top1new - twopi
      else if (top1new .le. -pi) then
         top1new = top1new + twopi
      end if

      top2new = - (thetanew - theta(neg_x(i)))
      if (top2new .gt. pi) then
         top2new = top2new-twopi
      else if (top2new .le. -pi) then
         top2new = top2new + twopi
      end if
      
      top3new = theta(pos_y(i)) - thetanew
      if (top3new .gt. pi) then
         top3new = top3new - twopi
      else if (top3new .le. -pi) then
         top3new = top3new + twopi
      end if

      top4new = - (theta(pos_x(i)) - thetanew)
      if (top4new .gt. pi) then
         top4new = top4new - twopi
      else if (top4new .le. -pi) then
         top4new = top4new + twopi
      end if

      ! METROPOLIS FILTER

      Uold = 0.5 * (top1old * top1old + top2old * top2old + top3old * top3old + top4old * top4old)
      Unew = 0.5 * (top1new * top1new + top2new * top2new + top3new * top3new + top4new * top4new)
      deltaU = Unew - Uold

      if ((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) then
         theta(i) = thetaNew
         top_x(i) = top1new
         top_y(i) = top2new
         top_x(pos_y(i)) = top3new
         top_y(pos_x(i)) = top4new
         accept = accept + 1
      end if
   end do
  
  return
end subroutine markov_chain_HXY
