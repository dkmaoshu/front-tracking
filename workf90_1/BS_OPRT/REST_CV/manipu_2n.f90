module manipulation_2
! This module is for fixing possible crossing of two stacked 
! discontinuity cells.

use discontinuity_curves

implicit none

public  manipulate_2

private memo_cleaning, truncation_1, distribution_1 


contains


subroutine manipulate_2 

implicit none

if(cvv(1)%status.eq.'asleep ') then
 print*, 'The curve list does not exist.'; return
end if

call memo_cleaning

call truncation_1

call distribution_1

end subroutine manipulate_2


subroutine memo_cleaning

implicit none

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
integer :: head_mark
logical :: head_switch
  
if(associated(cvv(1)%begin%next)) then
 temp=>cvv(1)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if
  
do while(associated(temp).and.head_switch) 

 g_cell=temp%g_cell
 temp%memo=0.0d0
 
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 endif
  
end do   

end subroutine memo_cleaning


subroutine truncation_1

implicit none

type(critical_cell), pointer :: temp, temps_l, temps_r
type(geo_info) :: g_cell, g_cells_l, g_cells_r
type(phy_info) :: p_cell, p_cells_l, p_cells_r
character*6 :: other_side
type(state) :: difference
logical :: left_truncate, right_truncate
integer :: head_mark
logical :: head_switch
  
if(associated(cvv(1)%begin%next)) then
 temp=>cvv(1)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if
  
do while(associated(temp).and.head_switch) 
 
 if(temp%l_stk.eq.adss_info(0,0).and.temp%r_stk.eq.adss_info(0,0)) then
  goto 10
 end if

 g_cell=temp%g_cell
 
 if(temp%l_stk.ne.adss_info(0,0)) then
  call visit(temp%l_stk,temps_l)
  call truncate(temps_l,'left  ',left_truncate,other_side)
 end if
 if(temp%r_stk.ne.adss_info(0,0)) then
  call visit(temp%r_stk,temps_r)
  call truncate(temps_r,'right ',right_truncate,other_side)
 end if
  
 if(left_truncate.and.right_truncate) call error_message
   
 10 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 endif
  
end do   


contains


subroutine truncate(temps,side,trunt,other_side)

implicit none

type(critical_cell), pointer :: temps
character*6, intent(in) :: side
logical, intent(out) :: trunt
character*6, intent(out) :: other_side

type(phy_info) :: p_cell, p_cells
real(8) :: ratio, ratios, difference
type(state) :: jump, jumps, aa, aas

trunt=.false.

if(temps%l_stk.eq.temp%address) then
 other_side='left  '
else
 if(temp%r_stk.eq.temp%address) then
  other_side='right '
 else
  call error_message
 end if
end if

p_cell=temp%p_cell; p_cells=temps%p_cell

select case(side) 
  
 case('left  ')
  jump=p_cell%l_state-p_cell%r_state
  aa=p_cell%or_state-p_cell%r_state

 case('right ')
  jump=p_cell%r_state-p_cell%l_state
  aa=p_cell%or_state-p_cell%l_state

 case default; call error_message

end select

ratio=aa%value/jump%value

select case(side) 
  
 case('left  ')
  jumps=p_cells%l_state-p_cells%r_state
  aas=p_cells%or_state-p_cells%r_state
   
 case('right ')
  jumps=p_cells%r_state-p_cells%l_state
  aas=p_cells%or_state-p_cells%l_state

 case default; call error_message

end select

ratios=aas%value/jumps%value

if(ratio+ratios.lt.1.0d0) then
 trunt=.true.
 difference=1.0-ratio-ratios
 temp%p_cell%or_state=temp%p_cell%or_state+difference*jump
 temp%memo=difference*jump
 temps%p_cell%or_state=temps%p_cell%or_state+difference*jumps
 temps%memo=difference*jumps
end if
 
end subroutine truncate


end subroutine truncation_1
   

subroutine distribution_1

implicit none

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(state) :: difference
integer :: head_mark
logical :: head_switch
  
if(associated(cvv(1)%begin%next)) then
 temp=>cvv(1)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if
  
do while(associated(temp).and.head_switch) 

 g_cell=temp%g_cell
 difference=temp%memo

 if(dabs(difference).lt.1.0d-7) goto 10

 if(associated(temp%previous)) then
  if(associated(temp%next)) then
   temp%previous%p_cell%or_state=temp%previous%p_cell%or_state- &
   0.5d0*difference
   temp%next%p_cell%or_state=temp%next%p_cell%or_state-0.5d0*difference
  else
   temp%previous%p_cell%or_state=temp%previous%p_cell%or_state-difference
  end if
 else
  if(associated(temp%next)) then
   temp%next%p_cell%or_state=temp%next%p_cell%or_state-difference
  else
   call error_message
  end if
 end if

 10 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 endif
  
end do   

end subroutine distribution_1


end module manipulation_2


