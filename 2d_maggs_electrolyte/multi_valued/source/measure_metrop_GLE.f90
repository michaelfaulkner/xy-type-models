! **************************************
! MEASURE STATE OF SYSTEM
! **************************************

subroutine measure
  use variables
  implicit none
  integer i, j, n
  real*8 potential, chargeNumber, storeEfield_x, storeEfield_y

  call measure_chargeDensity
  call measure_Esum

  chargeNumber = 0.0
  potential = 0.0

  do i = 0, sites - 1
     
     if (rho(i) .ne. 0) then
        chargeNumber = chargeNumber + 1.0
     end if
     storeEfield_x = Efield_x(i)
     storeEfield_y = Efield_y(i)
     potential = potential + 0.5 * (storeEfield_x * storeEfield_x + storeEfield_y * storeEfield_y)
     
  end do

  chargeNumber = chargeNumber / volume
  
  ! STORE SAMPLES DRAWN FROM MARKOV CHAIN
  
  write(10,100) Esum_x
  write(11,100) Esum_y
  write(12,100) chargeNumber
  write(13,100) potential

100 format(F16.8)

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
  do i = 0, sites - 1
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
  real*8 charge
  do i = 0, sites - 1
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
  character(100) filename
  real*8 chargeAcceptanceRate, auxFieldAcceptanceRate, TSFAcceptanceRate

  chargeAcceptanceRate = float(accept_charge) / (2 * measurements * volume)
  auxFieldAcceptanceRate = float(accept_aux_field) / (measurements * volume)
  TSFAcceptanceRate = float(accept_TSF) / measurements

  write (filename, '( "temp_eq_", F4.2,"//acceptance_rates_GLE_", I3.3, "x", I3.3, "_temp", F4.2, ".dat" )' ) &
       T, side, side, T
  open(unit = 300, file = filename)
  
  write(300, 100) chargeAcceptanceRate
  write(300, 100) auxFieldAcceptanceRate
  if (globalTSFon .eq. 1) then
     write(300, 100) TSFAcceptanceRate
  end if
  close(300)

100 format(F16.8)
  
  return
end subroutine output_acceptance_rates
