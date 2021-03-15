! **************************************
! SETS NEIGHBOURS WITH PERIODIC BOUNDARY CONDITIONS
! **************************************

subroutine PBC
  use variables
  implicit none
  integer i
  ! mod(i + side + 1,side) RATHER THAN mod(i + 1,side), ETC. BELOW AS mod(x,side)
  ! DOESN'T RETURN VALUES IN THE INTERVAL [0,side) FOR NEGATIVE x
  do i = 0,sites - 1
     pos_x(i) = i + mod(i + side + 1,side) - mod(i + side,side)
     neg_x(i) = i + mod(i + side - 1,side) - mod(i + side,side)
     pos_y(i) = i + (mod(int(i / side) + side + 1,side) - mod(int(i / side) + side,side)) * side
     neg_y(i) = i + (mod(int(i / side) + side - 1,side) - mod(int(i / side) + side,side)) * side
  end do
  return
end subroutine PBC

! **************************************
! INITIAL SPIN CONFIGURATION
! **************************************

subroutine initial_spins(start)
  use variables
  implicit none
  integer start,i,j
  if (start.eq.0) then
     do i=0,sites-1
        theta(i) = 0.0
     end do
  else
     do i=0,sites-1
        theta(i) = twopi * rand()
     end do
  end if
  return
end subroutine initial_spins

! **************************************
! SET CHAIN LENGTH AND Nevents TO 0; CREATE NEW DIRECTORY AND FILES FOR NEW TEMP
! **************************************

subroutine initial_measure
use variables
implicit none
character(100) filename
character(100) temperature_directory
character(100), parameter :: temperature_string="/temp_eq_"

! OPENS NEW DIRECTORY IN WHICH TO SAVE THE MARKOV CHAIN FOR THE CURRENT TEMPERATURE
write (temperature_directory, '(A, F4.2)' ) trim(output_directory)//trim(temperature_string), T
call system ( 'mkdir -p ' // temperature_directory )

write (filename, '(A, F4.2, "//magn_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=10, file=filename)
write (filename, '(A, F4.2, "//magn_x_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=11,file=filename)
write (filename, '(A, F4.2, "//magn_y_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=12,file=filename)
write (filename, '(A, F4.2, "//cos_top_x_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=13,file=filename)
write (filename, '(A, F4.2, "//cos_top_y_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=14,file=filename)
write (filename, '(A, F4.2, "//sin_top_x_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=15,file=filename)
write (filename, '(A, F4.2, "//sin_top_y_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=16,file=filename)
write (filename, '(A, F4.2, "//Ebar_x_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=17,file=filename)
write (filename, '(A, F4.2, "//Ebar_y_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=18,file=filename)
write (filename, '(A, F4.2, "//vort_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=19,file=filename)
write (filename, '(A, F4.2, "//second_deriv_potential_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit = 20, file = filename)
write (filename, '(A, F4.2, "//potential_sample.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit = 21, file = filename)

return
end subroutine initial_measure