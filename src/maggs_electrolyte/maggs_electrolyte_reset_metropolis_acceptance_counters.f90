subroutine reset_metropolis_acceptance_counters
use variables
implicit none

no_of_accepted_field_rotations = 0
no_of_accepted_charge_hops = 0
no_of_accepted_external_global_moves = 0

end subroutine reset_metropolis_acceptance_counters
