program elementary_charge_maggs_electrolyte_algorithm
use variables
implicit none
character(100) :: config_file
integer i, j, k, seed, start
double precision Tincr

! verify that the something has been parsed to the exectuable
if (command_argument_count() /= 1) then
    write(6, *) 'Error: parse configuration file to executable'
    stop
end if
! read in config file
call get_command_argument(1, config_file)
open (unit=1, file=config_file)

call input(seed, start)
call PBC
call randinit(seed)
write(6, *) rand(seed)
  
T = Tmin
if (Tsteps == 0) then
    Tincr = 0.0
else
    Tincr = (Tmax - Tmin) / Tsteps
end if
  
do i = 0, Tsteps

    write(6, *) T
    beta = 1 / T
    if (T .eq. Tmin) then
        call initial_Efield
        call measure_chargeDensity
    end if
     
    do j = 1, thermSweeps
        call markov_chain_aux_field_GLE
        do k = 1, ratio_charge_updates
            call markov_chain_charges_GLE
        end do
        if (globalTSFon .eq. 1) then
            do k = 1, ratio_TSF_updates
                call markov_chain_TSF_GLE
            end do
        end if
    end do

    call initial_measure
     
    do j = 1, measurements
        call markov_chain_aux_field_GLE
        do k = 1, ratio_charge_updates
            call markov_chain_charges_GLE
        end do
        if (globalTSFon .eq. 1) then
            do k = 1, ratio_TSF_updates
                call markov_chain_TSF_GLE
            end do
        end if
        call measure
    end do

    call output_acceptance_rates
    T = T + Tincr
    width_of_proposal_interval = width_of_proposal_interval + magnitude_of_proposal_interval_increments

end do
end program elementary_charge_maggs_electrolyte_algorithm


! **************************************
! METROPOLIS MARKOV-CHAIN SUBROUTINE FOR CHARGES
! **************************************

subroutine markov_chain_charges_GLE
  use variables
  implicit none
  integer n, i, plusMinus, rhoInew, rhoIposNew
  double precision deltaU, EfieldOld, EfieldNew

   do n = 1, 2 * sites

      i = int(rand() * sites) + 1

      if (floor(2 * rand()) .eq. 0) then

         EfieldOld = Efield_x(i)
         plusMinus = 2 * floor(2 * rand()) - 1
         EfieldNew = EfieldOld + plusMinus * elementaryCharge
         deltaU = 0.5 * (EfieldNew * EfieldNew - EfieldOld * EfieldOld)
         rhoInew = rho(i) + plusMinus
         rhoIposNew = rho(pos_x(i)) - plusMinus
         
         ! METROPOLIS FILTER

         if (((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) .and. &
              ((rhoInew .eq. 0) .or. (abs(rhoInew) .eq. 1)) .and. &
                 ((rhoIposNew .eq. 0) .or. (abs(rhoIposNew) .eq. 1))) then
            Efield_x(i) = EfieldNew
            rho(i) = rhoINew
            rho(pos_x(i)) = rhoIposNew
            accept_charge = accept_charge + 1
         end if
            
      else

         EfieldOld = Efield_y(i)
         plusMinus = 2 * floor(2 * rand()) - 1
         EfieldNew = EfieldOld + plusMinus * elementaryCharge
         deltaU = 0.5 * (EfieldNew * EfieldNew - EfieldOld * EfieldOld)
         rhoInew = rho(i) + plusMinus
         rhoIposNew = rho(pos_y(i)) - plusMinus
         
         ! METROPOLIS FILTER

         if (((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) .and. &
              ((rhoInew .eq. 0) .or. (abs(rhoInew) .eq. 1)) .and. &
                 ((rhoIposNew .eq. 0) .or. (abs(rhoIposNew) .eq. 1))) then
            Efield_y(i) = EfieldNew
            rho(i) = rhoINew
            rho(pos_y(i)) = rhoIposNew
            accept_charge = accept_charge + 1
         end if

      end if

   end do
   return

 end subroutine markov_chain_charges_GLE

! **************************************
! METROPOLIS MARKOV-CHAIN SUBROUTINE FOR AUXILIARY FIELD
! **************************************

subroutine markov_chain_aux_field_GLE
  use variables
  implicit none
  integer n, i
  double precision thetaOld, thetaNew, deltaTheta, Uold, Unew, deltaU
  double precision Efield1old, Efield2old, Efield3old, Efield4old
  double precision Efield1new, Efield2new, Efield3new, Efield4new

   do n = 1, sites

      i = int(rand() * sites) + 1
      deltaTheta = width_of_proposal_interval * (rand() - 0.5)

      ! CALL OLD ELECTRIC FIELD

      Efield1old = Efield_x(i)
      Efield2old = Efield_y(i)
      Efield3old = Efield_x(pos_y(i))
      Efield4old = Efield_y(pos_x(i))

      ! PROPOSED ELECTRIC FIELD
      
      Efield1new = Efield1old + deltaTheta
      Efield2new = Efield2old - deltaTheta
      Efield3new = Efield3old - deltaTheta
      Efield4new = Efield4old + deltaTheta

      ! METROPOLIS FILTER

      Uold = 0.5 * (Efield1old * Efield1old + Efield2old * Efield2old + Efield3old * Efield3old + Efield4old * Efield4old)
      Unew = 0.5 * (Efield1new * Efield1new + Efield2new * Efield2new + Efield3new * Efield3new + Efield4new * Efield4new)
      deltaU = Unew - Uold

      if ((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) then
         Efield_x(i) = Efield1new
         Efield_y(i) = Efield2new
         Efield_x(pos_y(i)) = Efield3new
         Efield_y(pos_x(i)) = Efield4new
         accept_aux_field = accept_aux_field + 1
      end if
   end do
  
  return
end subroutine markov_chain_aux_field_GLE
