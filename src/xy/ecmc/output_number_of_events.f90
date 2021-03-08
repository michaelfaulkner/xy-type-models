! **************************************
! OUTPUT Nevents
! **************************************

subroutine output_Nevents
  use variables
  implicit none
  character(100) filename

  write (filename, '( "temp_eq_", F4.2,"//Nevents_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' )  T,side,side,T
  open(unit=20,file=filename)
  write(20,100) Nevents
  if (twist.eq.1) then
     write(20,100) accept_twist
  end if
  close(20)

100 format(I20)

  return
end subroutine output_Nevents
