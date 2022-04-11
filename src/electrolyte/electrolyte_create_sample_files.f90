subroutine create_sample_files(temperature_index)
use variables
implicit none
character(100) :: temperature_directory, filename
integer :: temperature_index, file_index

! opens new directory in which to save the sample files at the current temperature
write(temperature_directory, '(A, "/temp_", I2.2)') trim(output_directory), temperature_index
call system('mkdir -p ' // temperature_directory)

if (measure_electric_field_sum) then
    write(filename, '(A, "/temp_", I2.2, "/electric_field_sum.csv")') trim(output_directory), temperature_index
    open(unit=30, file=filename)
    call print_file_header(30)
end if

if (measure_potential) then
    write(filename, '(A, "/temp_", I2.2, "/potential.csv")') trim(output_directory), temperature_index
    open(unit=40, file=filename)
    call print_file_header(40)
end if

if (measure_external_global_moves) then
    write(filename, '(A, "/temp_", I2.2, "/external_global_moves.csv")') trim(output_directory), temperature_index
    open(unit=50, file=filename)
    call print_file_header(50)
end if

return
end subroutine create_sample_files


subroutine print_file_header(file_index)
use variables
implicit none
integer :: temperature_index, file_index

write(file_index, '("#", A30, "; ", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') &
        trim(algorithm_name), integer_lattice_length, integer_lattice_length, temperature

return
end subroutine print_file_header
