subroutine output_metropolis_acceptance_rates
use variables
implicit none
character(100), parameter :: temperature_string="/temp_eq_"
character(100) :: filename
double precision :: acceptance_rate_of_external_local_moves, acceptance_rate_of_external_global_moves

acceptance_rate_of_external_local_moves = float(no_of_accepted_local_moves) / (measurements * volume)
acceptance_rate_of_external_global_moves = float(no_of_accepted_external_global_moves) / (measurements * volume)

write (filename, '(A, F4.2, "//acceptance_rates.dat")' ) trim(output_directory)//trim(temperature_string), temperature
open(unit = 300, file = filename)

write(300, 100) acceptance_rate_of_external_local_moves
if (twist == 1) then
    write(300, 100) acceptance_rate_of_external_global_moves
end if
close(300)

100 format(F16.8)

return
end subroutine output_metropolis_acceptance_rates
