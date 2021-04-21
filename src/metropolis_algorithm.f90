program metropolis_algorithm
use variables
implicit none
integer :: i, j

call pre_simulation_processes

do i = 0, no_of_temperature_increments
    write(6, '(A, F4.2)') 'Temperature = ', temperature
    beta = 1.0 / temperature
    call create_sample_file

    do j = 1, therm_sweeps
        call metropolis_sweep
        if (twist == 1) then
            call attempt_external_global_move
        end if
        call draw_observations
    end do

    call reset_metropolis_acceptance_counters

    do j = 1, measurements
        call metropolis_sweep
        if (twist == 1) then
            call attempt_external_global_move
        end if
        call draw_observations
    end do

    call output_metropolis_acceptance_rates
    temperature = temperature + magnitude_of_temperature_increments
    width_of_proposal_interval = width_of_proposal_interval + magnitude_of_proposal_interval_increments
end do
end program metropolis_algorithm
