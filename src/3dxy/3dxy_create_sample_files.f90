subroutine create_sample_files(temperature_index)
use variables
implicit none
character(100) :: filename
integer :: temperature_index, file_index

if (measure_magnetisation) then
    write(filename, '(A, "/temp_", I2.2, "/magnetisation.csv")') trim(output_directory), temperature_index
    if (start_from_checkpoint) then
        call trim_existing_sample_file(filename, 30, temperature_index)
    else
        open(unit=30, file=filename)
        call print_file_header(30)
    end if
end if

if (measure_potential) then
    write(filename, '(A, "/temp_", I2.2, "/potential.csv")') trim(output_directory), temperature_index
    if (start_from_checkpoint) then
        call trim_existing_sample_file(filename, 40, temperature_index)
    else
        open(unit=40, file=filename)
        call print_file_header(40)
    end if
end if

return
end subroutine create_sample_files

subroutine print_file_header(file_index)
use variables
implicit none
integer :: temperature_index, file_index

write(file_index, '("#", A30, "; ", I0.4, "x", I0.4, "x", I0.4, " lattice sites; temperature = ", ES8.2)') &
        trim(algorithm_name), integer_lattice_length, integer_lattice_length, integer_lattice_length, temperature

return
end subroutine print_file_header
