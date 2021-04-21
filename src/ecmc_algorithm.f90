program ecmc_algorithm
use variables
implicit none
integer :: i, j

call pre_simulation_processes

do i = 0, no_of_temperature_increments
    write(6, '(A, F4.2)') 'Temperature = ', temperature
    beta = 1.0 / temperature
    call create_sample_file

    do j = 1, therm_sweeps
        call single_event_chain
        if (twist == 1) then
            call attempt_external_global_move
        end if
        call draw_observations
    end do

    no_of_events = 0
    no_of_accepted_external_global_moves = 0

    do j = 1, measurements
        call single_event_chain
        if (twist == 1) then
            call attempt_external_global_move
        end if
        call draw_observations
    end do

    call output_no_of_events
    temperature = temperature + magnitude_of_temperature_increments
end do
end program ecmc_algorithm


subroutine output_no_of_events
use variables
implicit none
character(100), parameter :: temperature_string="/temp_eq_"
character(100) :: filename

write (filename, '(A, F4.2, "//number_of_events.csv")' ) trim(output_directory)//trim(temperature_string), temperature
open(unit=300, file = filename)
if (twist /= 1) then
    write(300, 100) no_of_events
else
    write(300, 200) no_of_events, float(no_of_accepted_external_global_moves) / (measurements * volume)
end if
close(300)

100 format(I20)
200 format(I20, ", ", ES24.14)

return
end subroutine output_no_of_events
