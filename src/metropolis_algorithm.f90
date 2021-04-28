program metropolis_algorithm
use variables
implicit none
integer :: i, j
double precision :: start_time, end_time, acceptance_rate_of_field_rotations

call pre_simulation_processes
call cpu_time(start_time)

do i = 0, no_of_temperature_increments
    write(6, '(A, F4.2)') 'Temperature = ', temperature
    beta = 1.0d0 / temperature
    call create_sample_file

    call reset_metropolis_acceptance_counters
    do j = 1, no_of_equilibration_sweeps
        call metropolis_sweep
        if (use_external_global_moves) then
            call attempt_external_global_move
        end if
        ! step-size adaptor
        if (mod(j, 100) == 0) then
            acceptance_rate_of_field_rotations = 1.0d-2 * dfloat(no_of_accepted_field_rotations) / dfloat(no_of_sites) &
                                                    / (1.0d0 - charge_hop_proportion)
            if (acceptance_rate_of_field_rotations > target_acceptance_rate_of_field_rotations + 5.0d-2) then
                width_of_proposal_interval = 1.1d0 * width_of_proposal_interval
            else if (acceptance_rate_of_field_rotations < target_acceptance_rate_of_field_rotations - 5.0d-2) then
                width_of_proposal_interval = 9.0d-1 * width_of_proposal_interval
            end if
            no_of_accepted_field_rotations = 0
        end if
        call get_and_print_observation
    end do

    call reset_metropolis_acceptance_counters
    do j = 1, no_of_observations
        call metropolis_sweep
        if (use_external_global_moves) then
            call attempt_external_global_move
        end if
        call get_and_print_observation
    end do

    call output_metropolis_acceptance_rates
    temperature = temperature + magnitude_of_temperature_increments
end do

call cpu_time(end_time)
write(6, '(A, ES9.3, A)') 'Markov process complete; total runtime of the simulation = ', end_time - start_time, &
                            ' seconds.'

end program metropolis_algorithm
