! **************************************
! MEASURE STATE OF SYSTEM
! **************************************

subroutine measure
use variables
implicit none
integer i
double precision potential

potential = 0.0
do i = 1, sites
    potential = potential + 0.5 * (Efield_x(i) * Efield_x(i) + Efield_y(i) * Efield_y(i))
end do

call measure_Esum

! STORE SAMPLES DRAWN FROM MARKOV CHAIN
write(10, 100) Esum_x, Esum_y
write(11, 200) potential

100 format(ES24.14, ", ", ES24.14)
200 format(ES24.14)

return
end subroutine measure

! **************************************
! MEASURE NUMBER OF CHARGES
! **************************************

subroutine measure_Esum
  use variables
  implicit none
  integer i

  Esum_x = 0.0
  Esum_y = 0.0
  do i = 1, sites
     Esum_x = Esum_x + Efield_x(i)
     Esum_y = Esum_y + Efield_y(i)
  end do

  return
end subroutine measure_Esum

! **************************************
! MEASURE NUMBER OF CHARGES
! **************************************

subroutine measure_chargeDensity
  use variables
  implicit none
  integer i,j
  double precision charge
  do i = 1, sites
     charge = ( Efield_x(i) + Efield_y(i) - Efield_x(neg_x(i)) - Efield_y(neg_y(i)) ) / elementaryCharge
     rho(i) = floor(charge + 0.5)
  end do
  return
end subroutine measure_chargeDensity


! **************************************
! OUTPUT acceptance rates
! **************************************

subroutine output_acceptance_rates
use variables
implicit none
character(100), parameter :: temperature_string="/temp_eq_"
character(100) filename
double precision chargeAcceptanceRate, auxFieldAcceptanceRate, TSFAcceptanceRate

chargeAcceptanceRate = float(accept_charge) / (2 * measurements * volume)
auxFieldAcceptanceRate = float(accept_aux_field) / (measurements * volume)
TSFAcceptanceRate = float(accept_TSF) / measurements

write (filename, '(A, F4.2, "//acceptance_rates.dat")' ) trim(output_directory)//trim(temperature_string), T
open(unit=300, file=filename)

write(300, 100) auxFieldAcceptanceRate
if (globalTSFon == 1) then
    write(300, 100) TSFAcceptanceRate
end if
write(300, 100) chargeAcceptanceRate
close(300)

100 format(F16.8)

return
end subroutine output_acceptance_rates