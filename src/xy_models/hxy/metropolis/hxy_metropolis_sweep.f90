subroutine metropolis_sweep
use variables
implicit none
integer :: n, i
double precision :: thetaOld, thetaNew, deltaTheta, Uold, Unew, deltaU, top1new, top2new, top3new, top4new

do n = 1, sites
    i = int(rand() * sites)  ! todo remove pick and replace (everywhere)
    thetaOld = theta(i)
    deltaTheta = 2.0d0 * proposalInterval * (rand() - 0.5d0)
    thetaNew = mod(thetaOld + deltaTheta, twopi)

    top1new = modulo(thetanew - theta(neg_y(i)) + pi, twopi) - pi
    top2new = modulo(- thetanew + theta(neg_x(i)) + pi, twopi) - pi
    top3new = modulo(theta(pos_y(i)) - thetanew + pi, twopi) - pi
    top4new = modulo(- theta(pos_x(i)) + thetanew + pi, twopi) - pi

    Uold = 0.5d0 * (top_x(i) * top_x(i) + top_y(i) * top_y(i) + &
                        top_x(pos_y(i)) * top_x(pos_y(i)) + top_y(pos_x(i)) * top_y(pos_x(i)))
    Unew = 0.5d0 * (top1new * top1new + top2new * top2new + top3new * top3new + top4new * top4new)
    deltaU = Unew - Uold

    if ((deltaU < 0.0d0) .or. (rand() < exp(- beta * deltaU))) then
        theta(i) = thetaNew
        top_x(i) = top1new
        top_y(i) = top2new
        top_x(pos_y(i)) = top3new
        top_y(pos_x(i)) = top4new
        no_of_accepted_local_moves = no_of_accepted_local_moves + 1
    end if
end do

return
end subroutine metropolis_sweep
