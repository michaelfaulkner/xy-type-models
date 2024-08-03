subroutine process_sample(sample_index)
use variables
implicit none
integer :: i, sample_index, raw_electric_field_zero_mode(2), get_topological_sector_component, topological_sector(2)
double precision :: potential

if (measure_electric_field_sum) then
    if (print_samples) then
        write(30, 100) - two_pi * dfloat(net_charge_displacement(1)), - two_pi * dfloat(net_charge_displacement(2))
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        ! raw_electric_field_zero_mode = \sum_r E(r) / two_pi = no_of_sites * Ebar / two_pi = - net_charge_displacement
        raw_electric_field_zero_mode(1) = - dfloat(net_charge_displacement(1))
        raw_electric_field_zero_mode(2) = - dfloat(net_charge_displacement(2))
        raw_electric_field_zero_mode_sum(1) = raw_electric_field_zero_mode_sum(1) + raw_electric_field_zero_mode(1)
        raw_electric_field_zero_mode_sum(2) = raw_electric_field_zero_mode_sum(2) + raw_electric_field_zero_mode(2)
        raw_electric_field_zero_mode_squared_sum = raw_electric_field_zero_mode_squared_sum + &
                                            raw_electric_field_zero_mode(1) ** 2 + raw_electric_field_zero_mode(2) ** 2
        raw_electric_field_zero_mode_quartic_sum = raw_electric_field_zero_mode_quartic_sum + &
                                    (raw_electric_field_zero_mode(1) ** 2 + raw_electric_field_zero_mode(2) ** 2) ** 2

        topological_sector(1) = get_topological_sector_component(raw_electric_field_zero_mode(1))
        topological_sector(2) = get_topological_sector_component(raw_electric_field_zero_mode(2))
        topological_sector_sum(1) = topological_sector_sum(1) + topological_sector(1)
        topological_sector_sum(2) = topological_sector_sum(2) + topological_sector(2)
        topological_sector_squared_sum = topological_sector_squared_sum + &
                                            topological_sector(1) ** 2 + topological_sector(2) ** 2
        topological_sector_quartic_sum = topological_sector_quartic_sum + &
                                            (topological_sector(1) ** 2 + topological_sector(2) ** 2) ** 2
    end if
end if

if (measure_potential) then
    potential = 0.0d0
    do i = 1, no_of_sites
        potential = potential + 0.5d0 * (electric_field(i, 1) ** 2 + electric_field(i, 2) ** 2)
    end do
    if (print_samples) then
        write(40, 200) potential
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        potential_sum = potential_sum + potential
        potential_squared_sum = potential_squared_sum + potential ** 2
    end if
end if

if (measure_external_global_moves) then
    write(50, 300) external_global_moves(1), external_global_moves(2)
end if

100 format(ES25.14, ",", ES25.14)
200 format(ES25.14)
300 format(I10, ",", I10)

return
end subroutine process_sample
