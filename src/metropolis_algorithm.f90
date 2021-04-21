program metropolis_algorithm
use variables
implicit none
integer :: i, j
double precision :: acceptance_rate_of_field_rotations

call pre_simulation_processes

do i = 0, no_of_temperature_increments
    write(6, '(A, F4.2)') 'Temperature = ', temperature
    beta = 1.0d0 / temperature
    call create_sample_file

    call reset_metropolis_acceptance_counters
    do j = 1, therm_sweeps
        call metropolis_sweep
        if (use_external_global_moves) then
            call attempt_external_global_move
        end if
        ! step-size adaptor
        if (mod(j, 100) == 0) then
            acceptance_rate_of_field_rotations = 1.0d-2 * float(no_of_accepted_field_rotations) / volume
            if (acceptance_rate_of_field_rotations > target_acceptance_rate_of_field_rotations + 5.0d-2) then
                width_of_proposal_interval = 1.1d0 * width_of_proposal_interval
            else if (acceptance_rate_of_field_rotations < target_acceptance_rate_of_field_rotations - 5.0d-2) then
                width_of_proposal_interval = 9.0d-1 * width_of_proposal_interval
            end if
            no_of_accepted_field_rotations = 0
        end if
        call draw_observations
    end do

    call reset_metropolis_acceptance_counters
    do j = 1, measurements
        call metropolis_sweep
        if (use_external_global_moves) then
            call attempt_external_global_move
        end if
        call draw_observations
    end do

    call output_metropolis_acceptance_rates
    temperature = temperature + magnitude_of_temperature_increments
end do
end program metropolis_algorithm
