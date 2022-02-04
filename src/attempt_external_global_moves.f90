subroutine attempt_external_global_moves
implicit none

if (floor(2.0d0 * rand()) == 0) then
    call attempt_single_external_global_move(1)
    call attempt_single_external_global_move(2)
else
    call attempt_single_external_global_move(2)
    call attempt_single_external_global_move(1)
end if

end subroutine attempt_external_global_moves
