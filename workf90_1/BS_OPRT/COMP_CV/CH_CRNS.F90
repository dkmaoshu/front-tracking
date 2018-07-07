module corner_check
! This module is for the corner-check of the 'xy'-type critcal
! cell.

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

implicit none

type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell, gp_cell, gn_cell

public check_corners
private ch_corners


contains


subroutine check_corners
! This subroutine checks 'xy'-type critical cells on auxiliary
! discontinuity curves to see if it is a corner of the curves

implicit none

integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake') then
  call ch_corners(acvv(i))
 end if
end do

end subroutine check_corners


subroutine ch_corners(acv)

implicit none

type(auxi_discv), intent(in) :: acv

integer :: head_mark
logical :: head_switch
character*3 :: corner, smooth

nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch)
 
 corner='no '
 g_cell=temp%g_cell
 if(g_cell%g_type.eq.'xy ') then
  if(.not.associated(temp%previous)) then
   corner='yes'
  else
   if(.not.associated(temp%next)) then
    corner='yes'
   else
    gp_cell=temp%previous%g_cell; gn_cell=temp%next%g_cell
    if(gp_cell%g_type.eq.'xy '.or.gn_cell%g_type.eq.'xy ') then
     corner='yes'
    else
     call check_smooth(gp_cell,gn_cell,smooth)
     if(smooth.eq.'no ') corner='yes'
    end if
   end if
  end if
 end if
 
 temp%corner=corner

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine ch_corners

end module corner_check 