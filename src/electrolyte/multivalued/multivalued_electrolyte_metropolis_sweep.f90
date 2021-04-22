subroutine metropolis_sweep
use variables
implicit none
integer i

call markov_chain_aux_field_GLE
do i = 1, ratio_charge_updates
    call markov_chain_charges_GLE
end do

return
end subroutine metropolis_sweep


subroutine markov_chain_charges_GLE
  use variables
  implicit none
  integer n, i, plusMinus
  double precision deltaU, EfieldOld, EfieldNew

   do n = 1, 2 * no_of_sites

      i = int(rand() * no_of_sites) + 1

      if (floor(2 * rand()) .eq. 0) then

         EfieldOld = electric_field_x(i)
         plusMinus = 2 * floor(2 * rand()) - 1
         EfieldNew = EfieldOld + plusMinus * elementary_charge
         deltaU = 0.5 * (EfieldNew * EfieldNew - EfieldOld * EfieldOld)
         
         ! METROPOLIS FILTER

         if ((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) then
            electric_field_x(i) = EfieldNew
            no_of_accepted_charge_hops = no_of_accepted_charge_hops + 1
         end if
            
      else

         EfieldOld = electric_field_y(i)
         plusMinus = 2 * floor(2 * rand()) - 1
         EfieldNew = EfieldOld + plusMinus * elementary_charge
         deltaU = 0.5 * (EfieldNew * EfieldNew - EfieldOld * EfieldOld)

         ! METROPOLIS FILTER

         if ((deltaU .lt. 0.0) .or. (exp(- beta * deltaU) .gt. rand())) then
            electric_field_y(i) = EfieldNew
            no_of_accepted_charge_hops = no_of_accepted_charge_hops + 1
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

   do n = 1, no_of_sites

      i = int(rand() * no_of_sites) + 1
      deltaTheta = width_of_proposal_interval * (rand() - 0.5)

      ! CALL OLD ELECTRIC FIELD

      Efield1old = electric_field_x(i)
      Efield2old = electric_field_y(i)
      Efield3old = electric_field_x(pos_y(i))
      Efield4old = electric_field_y(pos_x(i))

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
         electric_field_x(i) = Efield1new
         electric_field_y(i) = Efield2new
         electric_field_x(pos_y(i)) = Efield3new
         electric_field_y(pos_x(i)) = Efield4new
         no_of_accepted_field_rotations = no_of_accepted_field_rotations + 1
      end if
   end do
  
  return
end subroutine markov_chain_aux_field_GLE
