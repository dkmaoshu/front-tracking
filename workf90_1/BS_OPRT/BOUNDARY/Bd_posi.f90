module compute_boundary_positions
! This module is for the implementation of boundary conditions for
! auxiliary discontinuity curves. Currently, only infinitive
! boundary condition is implemented.

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

implicit none
type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell

public boundary_positions
private temp, g_cell, compute_bdp


contains


subroutine boundary_positions
! Implement all.
implicit none
integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake'.and.acvv(i)% &
    cv_type.ne.'circular') then
  call compute_bdp(acvv(i))
 end if
end do

end subroutine boundary_positions


subroutine compute_bdp(acv)
! implement a single.

implicit none
type(auxi_discv), intent(inout) :: acv

if(acv%begin_end.eq.0.and.associated(acv%begin%next)) then
 temp=>acv%begin%next
 g_cell=temp%g_cell
 g_cell%dis_posi(1)=g_cell%dis_posi(2)
 temp%g_cell=g_cell
end if

if(acv%end_end.eq.0.and.associated(acv%eend%previous)) then
 temp=>acv%eend%previous
 g_cell=temp%g_cell
 g_cell%dis_posi(2)=g_cell%dis_posi(1)
 temp%g_cell=g_cell
end if

end subroutine compute_bdp

end module compute_boundary_positions