subroutine process_sample(sample_index)
use variables
implicit none
integer :: i, sample_index, get_topological_sector_component, topological_sector(2)
double precision :: potential, raw_electric_field_zero_mode(2)

if (measure_electric_field_sum) then
    if (print_samples) then
        write(30, 100) - two_pi * dfloat(net_charge_displacement(1)), - two_pi * dfloat(net_charge_displacement(2))
    end if
    if (sample_index >= no_of_equilibration_sweeps) then
        ! raw_electric_field_zero_mode = \sum_r E(r) / no_of_sites / two_pi = Ebar / two_pi
        !                                                                   = - net_charge_displacement / no_of_sites
        ! nb, we divide by no_of_sites at this point to avoid infinities at large system size and long timescale
        raw_electric_field_zero_mode(1) = - dfloat(net_charge_displacement(1)) / dfloat(no_of_sites)
        raw_electric_field_zero_mode(2) = - dfloat(net_charge_displacement(2)) / dfloat(no_of_sites)
        raw_electric_field_zero_mode_sum(1) = raw_electric_field_zero_mode_sum(1) + raw_electric_field_zero_mode(1)
        raw_electric_field_zero_mode_sum(2) = raw_electric_field_zero_mode_sum(2) + raw_electric_field_zero_mode(2)
        raw_electric_field_zero_mode_squared_sum = raw_electric_field_zero_mode_squared_sum + &
                                    raw_electric_field_zero_mode(1) ** 2 + raw_electric_field_zero_mode(2) ** 2
        raw_electric_field_zero_mode_quartic_sum = raw_electric_field_zero_mode_quartic_sum + &
                                    (raw_electric_field_zero_mode(1) ** 2 + raw_electric_field_zero_mode(2) ** 2) ** 2

        topological_sector(1) = get_topological_sector_component(net_charge_displacement(1))
        topological_sector(2) = get_topological_sector_component(net_charge_displacement(2))
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
        potential_quartic_sum = potential_quartic_sum + potential ** 4
    end if
end if

if ((measure_external_global_moves).and.(print_samples)) then
    write(50, 300) external_global_moves(1), external_global_moves(2)
end if

100 format(ES25.14, ",", ES25.14)
200 format(ES25.14)
300 format(I10, ",", I10)

return
end subroutine process_sample


function get_topological_sector_component(net_charge_displacement_component)
use variables
implicit none
integer :: get_topological_sector_component, net_charge_displacement_component
! *** NB, this differs from the XY models version as it takes an integer ***
! w_{x / y} = floor((sum_r E_{x / y}(r) + pi L) / (2pi L)), where w \in Z^2 is the topological_sector
! to compute this accurately w/small epsilon > 0: w_{x / y} = floor((sum_r E_{x / y}(r) + pi L + 2pi epsilon) / (2pi L))
! this ensures edge cases are correctly processed, eg, Ebar_x = - pi / L does not incorrectly set w_x = -1
! *** NB, the epsilon was mistakenly absent in the code of Phys. Rev. B 91, 155412 (2015) ***
! *** this led to lower quality chi_w estimates ***
get_topological_sector_component = floor(-dfloat(net_charge_displacement_component) / dfloat(integer_lattice_length) + &
                                            0.5d0 + 1.0d-8 / dfloat(integer_lattice_length))
end function get_topological_sector_component
