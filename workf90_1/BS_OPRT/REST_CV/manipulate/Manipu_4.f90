module manipulation_4
! After the reset of discontinuity curve some stack relation may be 
! missing. This module is for fixing this problem. This fixing is 
! tentatively, a through solution will introducing four-neighbor stack
! relation.

use auxiliary_discontinuity_curves

implicit none

public  manipulate_4
private fix_curve


contains


subroutine manipulate_4 

implicit none

if(acvv(1)%status.eq.'asleep ') then
 print*, 'The curve list does not exist.'; return
end if

call fix_curve

end subroutine manipulate_4


subroutine fix_curve

implicit none

type(auxi_crit_cell), pointer :: temp, temp_s, boundary_end
integer :: head_mark, end_mark
logical :: head_switch

if(associated(acvv(1)%begin%next)) then
 temp=>acvv(1)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
 boundary_end=>acvv(1)%eend_boundary%previous
 if(acvv(1)%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acvv(1)%eend%previous%next%address%idx
 end if
else
 print*, 'The curve list is empty.'; return
end if
  
do while(associated(temp).and.head_switch) 

 if(temp%l_stk.ne.adss_info(0,0)) then
  call visit(temp%l_stk,temp_s)
  if(temp_s%l_stk.ne.temp%address.and.temp_s%r_stk.ne.temp%address) then
   call fix_stack('left  ')
  end if
 end if

 if(temp%r_stk.ne.adss_info(0,0)) then
  call visit(temp%r_stk,temp_s)
  if(temp_s%l_stk.ne.temp%address.and.temp_s%r_stk.ne.temp%address) then
   call fix_stack('right ')
  end if
 end if

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 endif
  
end do 


contains  


subroutine fix_stack(side)

implicit none

character*6, intent(in) :: side

type(geo_info) :: g_cell, g_cell_s
integer :: corner
character*6 :: other_side

g_cell=temp%g_cell; g_cell_s=temp_s%g_cell
other_side=side_shift(side)
do corner=1, 4
 if(g_cell%point(corner).eq.other_side) exit
end do
if(corner.gt.4) call error_message
select case(g_cell_s%point(corner))
 case('left  '); temp_s%l_stk=temp%address
 case('right '); temp_s%r_stk=temp%address
 case default; call error_message
end select

end subroutine fix_stack


end subroutine fix_curve


end module manipulation_4


