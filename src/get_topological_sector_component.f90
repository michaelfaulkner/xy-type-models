function get_topological_sector_component(raw_electric_field_zero_mode_component)
use variables
implicit none
integer :: get_topological_sector_component
double precision :: raw_electric_field_zero_mode_component
! w_{x / y} = floor((sum_r E_{x / y}(r) + pi L) / (2pi L)), where w \in Z^2 is the topological_sector
! compute this accurately w/small epsilon > 0: w_{x / y} = floor((sum_r E_{x / y}(r) + pi L - 2pi epsilon) / (2pi L))
! get_topological_sector_component = floor((sum_of_emergent_field_component + pi * dfloat(integer_lattice_length) - &
!                                             two_pi * 1.0d-8) / (two_pi * dfloat(integer_lattice_length)))
get_topological_sector_component = floor(raw_electric_field_zero_mode_component / two_pi / &
                                    dfloat(integer_lattice_length) + 0.5d0 - 1.0d-8 / dfloat(integer_lattice_length))
end function get_topological_sector_component
