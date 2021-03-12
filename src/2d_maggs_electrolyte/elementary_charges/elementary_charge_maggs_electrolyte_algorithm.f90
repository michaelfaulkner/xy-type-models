! **************************************
! Metropolis sampling from GLE Coulomb distribution (elementary charges)
! **************************************

! **************************************
! Main program
! **************************************

program metrop_GLE

  use variables
  implicit none
  integer i, j, k, seed, start
  real*8 Tincr

  open (unit = 1, file = 'initial.in')
  call input(seed, start)
  call PBC
  call randinit(seed)
  write(6,*) rand(seed)
  
  T = Tmin
  if (Tsteps .eq. 0) then
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
     
     do j = 0, thermSweeps - 1
        call markov_chain_aux_field_GLE
        do k = 0, ratio_charge_updates - 1
           call markov_chain_charges_GLE
        end do
        if (globalTSFon .eq. 1) then
           do k = 0, ratio_TSF_updates - 1
              call markov_chain_TSF_GLE
           end do
        end if
     end do

     call initial_measure
     
     do j = 0, measurements - 1
        
        call markov_chain_aux_field_GLE
        do k = 0, ratio_charge_updates - 1
           call markov_chain_charges_GLE
        end do
        if (globalTSFon .eq. 1) then
           do k = 0, ratio_TSF_updates - 1
              call markov_chain_TSF_GLE
           end do
        end if

        call measure

     end do

     call output_acceptance_rates

     T = T + Tincr
     proposalInterval = proposalInterval + deltaProposalInterval
  end do
end program metrop_GLE


! **************************************
! METROPOLIS MARKOV-CHAIN SUBROUTINE FOR CHARGES
! **************************************

subroutine markov_chain_charges_GLE
  use variables
  implicit none
  integer n, i, plusMinus, rhoInew, rhoIposNew
  real*8 deltaU, EfieldOld, EfieldNew

   do n = 0, 2 * sites - 1

      i = int(rand() * sites)

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
  real*8 thetaOld, thetaNew, deltaTheta, Uold, Unew, deltaU
  real*8 Efield1old, Efield2old, Efield3old, Efield4old
  real*8 Efield1new, Efield2new, Efield3new, Efield4new

   do n = 0, sites - 1

      i = int(rand() * sites)
      deltaTheta = 2. * proposalInterval * (rand() - 0.5)

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
