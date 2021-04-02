subroutine metropolis_sweep
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
         no_of_accepted_local_moves = no_of_accepted_local_moves + 1
      end if
   end do
  
  return
end subroutine metropolis_sweep
