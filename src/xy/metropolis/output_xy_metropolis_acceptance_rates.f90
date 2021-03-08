! **************************************
! OUTPUT acceptance rates
! **************************************

subroutine output_acceptance_rates
  use variables
  implicit none
  character(100) filename
  real*8 acceptanceRate, twistAcceptanceRate

  acceptanceRate = float(accept) / (measurements * volume)
  twistAcceptanceRate = float(accept_twist) / (measurements * volume)

  write (filename, '( "temp_eq_", F4.2,"//acceptance_rates_XY_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' ) &
       T, side, side, T
  open(unit = 300, file = filename)

  write(300, 100) acceptanceRate
  if (twist .eq. 1) then
     write(300, 100) twistAcceptanceRate
  end if
  close(300)

100 format(F16.8)

  return
end subroutine output_acceptance_rates
