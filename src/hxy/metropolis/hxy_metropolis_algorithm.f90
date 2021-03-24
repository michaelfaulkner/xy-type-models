program hxy_metropolis_algorithm
use variables
implicit none
character(100) :: config_file
integer :: i, j, seed, start
double precision :: Tincr

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
    beta = 1.0d0 / T

    do j = 0, thermSweeps - 1
        call markov_chain_HXY
    end do

    accept = 0
    accept_twist = 0
    call initial_measure

    do j = 0, measurements - 1
        call markov_chain_HXY
        if (twist == 1) then
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
integer :: n, i
double precision :: thetaOld, thetaNew, deltaTheta, Uold, Unew, deltaU, top1new, top2new, top3new, top4new

do n = 1, sites
    i = int(rand() * sites)  ! todo remove pick and replace (everywhere)
    thetaOld = theta(i)
    deltaTheta = 2.0d0 * proposalInterval * (rand() - 0.5d0)
    thetaNew = modulo(thetaOld + deltaTheta + pi, twopi) - pi

    top1new = modulo(thetanew - theta(neg_y(i)) + pi, twopi) - pi
    top2new = modulo(- thetanew + theta(neg_x(i)) + pi, twopi) - pi
    top3new = modulo(theta(pos_y(i)) - thetanew + pi, twopi) - pi
    top4new = modulo(- theta(pos_x(i)) + thetanew + pi, twopi) - pi

    Uold = 0.5d0 * (top_x(i) * top_x(i) + top_y(i) * top_y(i) + &
                        top_x(pos_y(i)) * top_x(pos_y(i)) + top_y(pos_x(i)) * top_y(pos_x(i)))
    Unew = 0.5d0 * (top1new * top1new + top2new * top2new + top3new * top3new + top4new * top4new)
    deltaU = Unew - Uold

    if ((deltaU < 0.0d0) .or. (rand() < exp(- beta * deltaU))) then
        theta(i) = modulo(thetaNew + pi, twopi) - pi
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
character(100) :: filename
double precision :: acceptanceRate, twistAcceptanceRate

acceptanceRate = float(accept) / (measurements * volume)
twistAcceptanceRate = float(accept_twist) / (measurements * volume)

write (filename, '(A, F4.2, "//acceptance_rates.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=300, file=filename)

write(300, 100) acceptanceRate
if (twist == 1) then
    write(300, 100) twistAcceptanceRate
end if
close(300)

100 format(F16.8)

return
end subroutine output_acceptance_rates
