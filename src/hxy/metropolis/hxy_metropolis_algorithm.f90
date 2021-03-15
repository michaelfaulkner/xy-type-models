program hxy_metropolis_algorithm
use variables
implicit none
character(100) :: config_file
integer i, j, seed, start
real*8 Tincr

! verify that the something has been parsed to the exectuable
if (command_argument_count() /= 1) then
    write(6, *) 'Error: parse configuration file to executable'
    stop
end if
! read in config file
call get_command_argument(1, config_file)
open (unit=1, file=config_file)

call input(seed, start)
call PBC
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
    beta = 1 / T
    if (T == Tmin) then
        call initial_spins(start)
    end if

    do j = 0, thermSweeps - 1
        call markov_chain_HXY
    end do

    accept = 0
    accept_twist = 0
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
end program hxy_metropolis_algorithm


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

! **************************************
! OUTPUT acceptance rates
! **************************************

subroutine output_acceptance_rates
use variables
implicit none
character(100), parameter :: temperature_string="/temp_eq_"
character(100) filename
real*8 acceptanceRate, twistAcceptanceRate

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
