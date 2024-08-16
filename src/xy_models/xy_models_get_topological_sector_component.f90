function get_topological_sector_component(non_normalised_emergent_field_zero_mode_component)
use variables
implicit none
integer :: get_topological_sector_component
double precision :: non_normalised_emergent_field_zero_mode_component
! *** NB, this differs from the electrolyte version as it takes a float ***
! w_{x / y} = floor((sum_r E_{x / y}(r) + pi L) / (2pi L)), where w \in Z^2 is the topological_sector
! compute this accurately w/small epsilon > 0: w_{x / y} = floor((sum_r E_{x / y}(r) + pi L + 2pi epsilon) / (2pi L))
! this ensures edge cases are correctly processed, eg, Ebar_x = - pi / L does not incorrectly set w_x = -1
! *** NB, the epsilon was mistakenly absent in the code of J. Phys.: Condens. Matter 29, 085402 (2017) ***
! *** this led to lower quality chi_w estimates ***
get_topological_sector_component = floor(non_normalised_emergent_field_zero_mode_component / two_pi / &
                                    dfloat(integer_lattice_length) + 0.5d0 + 1.0d-8 / dfloat(integer_lattice_length))
end function get_topological_sector_component
