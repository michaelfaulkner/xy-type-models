subroutine set_running_sums
use variables
implicit none
character(100) :: filename

if (measure_electric_field_sum) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/electric_field_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 200) raw_electric_field_zero_mode_sum(1)
        read(11, 200) raw_electric_field_zero_mode_sum(2)
        read(11, 200) raw_electric_field_zero_mode_squared_sum
        read(11, 200) raw_electric_field_zero_mode_quartic_sum
        read(11, 100) topological_sector_sum(1)
        read(11, 100) topological_sector_sum(2)
        read(11, 100) topological_sector_squared_sum
        read(11, 100) topological_sector_quartic_sum
        close(11)
    else
        raw_electric_field_zero_mode_sum = (/ 0.0d0, 0.0d0 /)
        raw_electric_field_zero_mode_squared_sum = 0.0d0
        raw_electric_field_zero_mode_quartic_sum = 0.0d0
        topological_sector_sum = (/ 0, 0 /)
        topological_sector_squared_sum = 0
        topological_sector_quartic_sum = 0
    end if
end if

if (measure_potential) then
    if (start_from_checkpoint) then
        write(filename, '(A, "/potential_running_sums.csv")') trim(output_directory)
        open(unit=11, file=filename)
        read(11, 200) potential_sum
        read(11, 200) potential_squared_sum
        read(11, 200) potential_quartic_sum
        close(11)
    else
        potential_sum = 0.0d0
        potential_squared_sum = 0.0d0
        potential_quartic_sum = 0.0d0
    end if
end if

100 format(I12)
200 format(ES25.14)

return
end subroutine set_running_sums
