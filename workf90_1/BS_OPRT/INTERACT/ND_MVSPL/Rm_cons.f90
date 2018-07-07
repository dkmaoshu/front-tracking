module corner_remove
! This module is written for removing corners.

use solution
! 'solution.f90'

implicit none

type(curve_plug), pointer :: tempp
integer :: nd_nmb, nd_nmb_n, cv_nmb
logical :: if_corner
character*6 :: other_end


public  remove_corners
private check_corner, tempp, cv_nmb, if_corner, nd_nmb, nd_nmb_n


contains


subroutine remove_corners

implicit none

integer :: head_mark
logical :: head_switch

! call check_ring(1,'counter  ')

do nd_nmb=1, ndn
 if(ndd(nd_nmb)%status.eq.'awake ') then
  tempp=>ndd(nd_nmb)%cv_rg%begin%next
  head_mark=tempp%address%idx; head_switch=.true.
  do while(head_switch)
   call check_corner
   if(if_corner) call remove_corner
   tempp=>tempp%next
   head_switch=(tempp%address%idx.ne.head_mark)
  end do
 end if
end do

end subroutine remove_corners


subroutine check_corner

implicit none

type(critical_cell), pointer :: temp
type(cv_plug_info) :: plug
type(geo_info) :: g_cell1, g_cell2, g_cell3
type(node_info) :: n_cell, n_cell_n
integer :: i3, j3

if_corner=.false.
n_cell=ndd(nd_nmb)%n_cell; other_end='      '; call clean_up_n_cell(n_cell_n)
plug=tempp%plug; cv_nmb=plug%cv_nmb
call clean_up_g_cell(g_cell1); call clean_up_g_cell(g_cell2); call clean_up_g_cell(g_cell3)

! Pick the first three critical cells.
select case(plug%end_type)

 case('begin ')
  temp=>cvv(cv_nmb)%begin
  if(associated(temp%next)) then
   temp=>temp%next; g_cell1=temp%g_cell
  else
   return
  end if
  if(associated(temp%next)) then
   temp=>temp%next; g_cell2=temp%g_cell
  else
   return
  end if
  if(associated(temp%next)) then
   temp=>temp%next; g_cell3=temp%g_cell; other_end='crit  '
  else
   if(cvv(cv_nmb)%end_end.le.0) call error_message
   nd_nmb_n=cvv(cv_nmb)%end_end; n_cell=ndd(nd_nmb_n)%n_cell
   other_end='node  '
  end if

 case('end   ')
  temp=>cvv(cv_nmb)%eend
  if(associated(temp%previous)) then
   temp=>temp%previous; g_cell1=temp%g_cell
  else
   return
  end if
  if(associated(temp%previous)) then
   temp=>temp%previous; g_cell2=temp%g_cell
  else
   return
  end if
  if(associated(temp%previous)) then
   temp=>temp%previous; g_cell3=temp%g_cell; other_end='crit  '
  else
   if(cvv(cv_nmb)%begin_end.le.0) call error_message
   nd_nmb_n=cvv(cv_nmb)%begin_end; n_cell=ndd(nd_nmb_n)%n_cell
   other_end='node  '
  end if

end select    

! Check whether there is a corner.
if(g_cell1%g_type.ne.'xy '.or.g_cell2%g_type.ne.'xy ') return
select case(other_end)
 case('crit  ')
  i3=g_cell3%x_idx; j3=g_cell3%y_idx
 case('node  ')
  i3=n_cell_n%x_idx; j3=n_cell_n%y_idx
end select 
if(i3.eq.n_cell%x_idx.and.j3.eq.n_cell%y_idx) call error_message
if(iabs(i3-n_cell%x_idx)+iabs(j3-n_cell%y_idx).eq.1) if_corner=.true.

end subroutine check_corner


subroutine remove_corner

implicit none

type(critical_cell), pointer :: temp
type(cv_plug_info) :: plug
type(geo_info) :: g_cell
type(state) :: difference1, difference2
character*6 :: side1, side2, position
integer :: single

call clean_up_plug(plug); call clean_up_g_cell(g_cell)
side1='      '; side2='      '
plug=tempp%plug

select case(plug%end_type)
 case('begin '); temp=>cvv(plug%cv_nmb)%begin%next; position='after '
 case('end   '); temp=>cvv(plug%cv_nmb)%eend%previous; position='before'
end select
g_cell=temp%g_cell
call find_single(g_cell,single)
side1=g_cell%point(single)  
call remove(temp,position,side1,difference1)

g_cell=temp%g_cell
call find_single(g_cell,single)
side2=g_cell%point(single)
if(side1.ne.side2) call error_message
call remove(temp,position, side2,difference2)

temp%p_cell%or_state=temp%p_cell%or_state+difference1+difference2


end subroutine remove_corner


end module corner_remove