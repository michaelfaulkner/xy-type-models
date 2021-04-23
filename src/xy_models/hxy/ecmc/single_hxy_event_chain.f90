subroutine single_event_chain
use variables
implicit none
integer :: i, active_spin_index, vetoeing_spin_index
integer, dimension (1:4) :: neighbouring_spin_indices
double precision :: uphill_distance_through_potential_space_before_next_event, shortest_distance_to_next_factor_event
double precision :: initial_two_spin_potential, final_two_spin_potential, distance_left_before_next_observation
double precision :: active_spin_value, non_active_spin_value, initial_spin_value_difference, final_spin_value_difference
double precision :: distance_to_next_factor_event, no_of_complete_spin_rotations

active_spin_index = int(dfloat(no_of_sites) * rand()) + 1
distance_left_before_next_observation = spin_space_distance_between_observations
! iterate until total distance covered in spin space reaches spin_space_distance_between_observations
do
    active_spin_value = spin_field(active_spin_index)
    neighbouring_spin_indices(1) = get_north_neighbour(active_spin_index)
    neighbouring_spin_indices(2) = get_south_neighbour(active_spin_index)
    neighbouring_spin_indices(3) = get_east_neighbour(active_spin_index)
    neighbouring_spin_indices(4) = get_west_neighbour(active_spin_index)

    shortest_distance_to_next_factor_event = 1.0d10
    ! iterate over neighbouring_spin_indices
    do i = 1, 4
        non_active_spin_value = spin_field(neighbouring_spin_indices(i))
        initial_spin_value_difference = modulo(active_spin_value - non_active_spin_value + pi, twopi) - pi
        uphill_distance_through_potential_space_before_next_event = - temperature * log(1.0d0 - rand())
        ! if factor derivative > 0
        if (initial_spin_value_difference > 0.0d0) then
            initial_two_spin_potential = 0.5d0 * initial_spin_value_difference * initial_spin_value_difference
            no_of_complete_spin_rotations = int((initial_two_spin_potential + &
                                                    uphill_distance_through_potential_space_before_next_event) / &
                                                pi_squared_over_two)
            final_two_spin_potential = (no_of_complete_spin_rotations + 1.0d0) * pi_squared_over_two - &
                                            initial_two_spin_potential - &
                                            uphill_distance_through_potential_space_before_next_event
            final_spin_value_difference = pi - sqrt(pi_squared - 2.0d0 * final_two_spin_potential)
            distance_to_next_factor_event = (no_of_complete_spin_rotations + 0.5d0) * twopi - &
                                                initial_spin_value_difference - final_spin_value_difference
        ! else: factor derivative < 0 ==> go to bottom of potential well
        else
            no_of_complete_spin_rotations = int(uphill_distance_through_potential_space_before_next_event / &
                                                    pi_squared_over_two)
            final_two_spin_potential = (no_of_complete_spin_rotations + 1.0d0) * pi_squared_over_two - &
                                            uphill_distance_through_potential_space_before_next_event
            final_spin_value_difference = pi - sqrt(pi_squared - 2.0d0 * final_two_spin_potential)
            distance_to_next_factor_event = (no_of_complete_spin_rotations + 0.5d0) * twopi - &
                                                initial_spin_value_difference - final_spin_value_difference
        end if

        if (distance_to_next_factor_event < shortest_distance_to_next_factor_event) then
            shortest_distance_to_next_factor_event = distance_to_next_factor_event
            vetoeing_spin_index = neighbouring_spin_indices(i)
        end if
    end do

    if (distance_left_before_next_observation < shortest_distance_to_next_factor_event) then
        ! update active spin value and exit event chain in order to observe the system
        spin_field(active_spin_index) = mod(active_spin_value + distance_left_before_next_observation, twopi)
        exit
    else
        ! update active spin value and continute event chain
        spin_field(active_spin_index) = mod(active_spin_value + shortest_distance_to_next_factor_event, twopi)
        distance_left_before_next_observation = distance_left_before_next_observation - &
                                                    shortest_distance_to_next_factor_event
        active_spin_index = vetoeing_spin_index
        no_of_events = no_of_events + 1
    end if
end do

return
end subroutine single_event_chain
