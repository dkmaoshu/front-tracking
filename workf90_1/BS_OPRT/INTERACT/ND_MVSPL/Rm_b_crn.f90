module remove_broken_corners
! This module is for remove of broken corners around an old node cell. The remove is
! implemented before treating the triple or multiple node cells.

use tools_4_connecting_n_triple
! 'tools_ct.f90'

implicit none

type(critical_cell), pointer :: temp, temp_second
character*3 :: if_corner, may_affect, if_empty
type(cv_plug_info) :: plug
type(geo_info) :: g_cell, g_cell_s
type(phy_info) :: p_cell_s
integer :: x_idx_s, y_idx_s 

public  remove_on_all_curves
private corner_check, affect_check, forced_remove, plug, g_cell, g_cell_s, x_idx_s, &
        y_idx_s, temp_second, temp, p_cell_s

contains


subroutine remove_on_all_curves

implicit none

type(curve_plug), pointer :: tempp_1
type(adss_plug) :: adss, adss_1 
					
tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)

 plug=tempp%plug; adss=tempp%address
! First, it should check if the disconnection is a corner and
! treat it accordingly.
! The disconnection is found to be a corner and the corresponding
! critical cell should be removed from the curve list.
 call corner_check

 if(if_corner.eq.'yes') then
! First, we should check if the following curves in pluging-ring
! have the same type of corner disconnnection, and the current
! remove will effect the following one. If it is the case, the 
! remove job should starts from the last curve of this kind in the 
! plugingring.

  call affect_check 
  select case(may_affect)
   case('no ')
    call forced_remove
   case('yes')
    head_mark=tempp%address%idx; head_switch=.true.
    tempp_1=>tempp%next 
    plug=tempp_1%plug; adss_1=tempp_1%address
    do while(head_switch)
     call corner_check
     if(if_corner.eq.'no ') then
      tempp_1=>tempp_1%previous; adss_1=tempp_1%address; exit
     end if
     call affect_check
     if(may_affect.eq.'no ') exit
     tempp_1=>tempp_1%next;
     adss_1=tempp_1%address; plug=tempp_1%plug
     head_switch=(adss_1%idx.ne.head_mark)
    end do
	tempp=>tempp_1

	plug=tempp_1%plug; adss=tempp_1%address
    do while(adss.ne.adss_1)
	 call corner_check
     call forced_remove
     tempp_1=>tempp_1%previous; adss=tempp_1%address
    end do
   end select
 end if

 tempp=>tempp%next
 head_switch=(tempp%address%idx.ne.head_mark)
 
end do

end subroutine remove_on_all_curves


subroutine forced_remove
! Remove the head critical cell and update the second critical cell
! in the list to fix the connection between the curve list and the 
! node cell.

implicit none

type(state) :: difference
character*3 :: direction
integer :: single, corner, the_end, edge_plus

! type(phy_info) :: ppp

! First, update the geometrical information in the second critical
! cell.
call find_single(g_cell,single)
select case(plug%end_type)
 case('begin '); the_end=1
 case('end   '); the_end=2
 case default; print*, 'Something is wrong!!!'; pause
end select
g_cell_s%edge(the_end)=g_cell%edge(the_end)
if(x_idx_s.ne.g_cell%x_idx) direction='xx '
if(y_idx_s.ne.g_cell%y_idx) direction='yy '
corner=point_shift(single,direction)
edge_plus=cycle_4(g_cell_s%edge(the_end)+1)
if(corner.ne.g_cell_s%edge(the_end).and.corner.ne.edge_plus) then
 print*, 'Something is wrong!!!'; pause
end if
select case(g_cell_s%edge(the_end))
 case(1,2)
  if(corner.eq.g_cell_s%edge(the_end)) then
   g_cell_s%dis_posi(the_end)=-0.4999999d0
  else
   g_cell_s%dis_posi(the_end)=0.4999999d0
  end if
 case(3,4)
  if(corner.eq.g_cell_s%edge(the_end)) then
   g_cell_s%dis_posi(the_end)=0.4999999d0
  else
   g_cell_s%dis_posi(the_end)=-0.4999999d0
  end if
end select
call find_type_from_edges(g_cell_s)
call find_points_from_edges(g_cell_s)
temp_second%g_cell=g_cell_s

! ppp=temp%p_cell

! Second, remove the first critical list, update the gird_map and
! compute the conservation difference caused.
difference=error_data
select case(plug%end_type)
 case('begin ')
  call remove(temp,'after ',g_cell%point(single),difference)
 case('end   ')
  call remove(temp,'before',g_cell%point(single),difference)
end select
ggd_cell(x_idx_s,y_idx_s)%ccpt(corner)%side=g_cell_s%point(corner)

! third, update the physical information in second critical cell.
p_cell_s%or_state=p_cell_s%or_state+difference
temp_second%p_cell=p_cell_s

! Finally, update the plug in the node cell and the grid map.
plug%edge=neighbor_edge(g_cell_s%edge(the_end))

! call check_list_c(3,'up    ')

end subroutine forced_remove


subroutine affect_check
! This subroutine checks whether the current remove may affect
! the following curve in the pluging-ring.

implicit none

type(curve_plug), pointer :: tempp_nxt
type(critical_cell), pointer :: temp_nxt
type(cv_plug_info) :: plug_nxt
type(geo_info) :: g_nxt, g_nxt_snd
integer :: x_dif, y_dif, x_dif_s, y_dif_s, max

may_affect='no '
x_dif=g_cell%x_idx-x_nd; y_dif=g_cell%y_idx-y_nd
x_dif_s=g_cell_s%x_idx-x_nd; y_dif_s=g_cell_s%y_idx-y_nd

nullify(temp_nxt); call clean_up_g_cell(g_nxt); call clean_up_g_cell(g_nxt_snd)
tempp_nxt=>tempp%next; plug_nxt=tempp_nxt%plug
call locate_head(plug_nxt,'first ',temp_nxt,if_empty)
if(associated(temp_nxt)) then
 g_nxt=temp_nxt%g_cell
 if(g_cell%x_idx.ne.g_nxt%x_idx.or.g_cell%y_idx.ne.g_nxt%y_idx) return
 nullify(temp_nxt)
 call locate_head(plug_nxt,'second',temp_nxt,if_empty)
 if(associated(temp_nxt)) then
  g_nxt_snd=temp_nxt%g_cell
  if(g_cell_s%x_idx.ne.g_nxt_snd%x_idx.or.g_cell_s%y_idx.ne.g_nxt_snd &
     %y_idx) return
 end if
else
 return
end if

if(x_dif*x_dif_s.le.0.and.y_dif*y_dif_s.le.0) then
 print*, 'Something is wrong!!!'; pause
end if
max=max0(iabs(x_dif),iabs(x_dif_s),iabs(y_dif),iabs(y_dif_s))
if(max.gt.1) then
 print*, 'Something is wrong!!!'; pause
end if

if(x_dif*x_dif_s.gt.0) then
 if(x_dif_s.gt.0) then
  if(y_dif.lt.0) may_affect='yes'
 else
  if(y_dif.gt.0) may_affect='yes'
 end if
else
 if(y_dif_s.gt.0) then
  if(x_dif.gt.0) may_affect='yes'
 else
  if(x_dif.lt.0) may_affect='yes'
 end if
end if

end subroutine affect_check


subroutine corner_check
! This subroutine checks if there is a corner disconnection.

implicit none

character*3 :: if_second_empty

nullify(temp_second); if_corner='no '
call clean_up_g_cell(g_cell_s)

call locate_head(plug,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
g_cell=temp%g_cell
if(iabs(x_nd-g_cell%x_idx)+iabs(y_nd-g_cell%y_idx).gt.1) then
 call locate_head(plug,'second',temp_second,if_second_empty)
 if(associated(temp_second)) then
  g_cell_s=temp_second%g_cell; p_cell_s=temp_second%p_cell
  x_idx_s=g_cell_s%x_idx; y_idx_s=g_cell_s%y_idx
  if(iabs(x_nd-x_idx_s)+iabs(y_nd-y_idx_s).le.1) if_corner='yes'
 end if
end if

end subroutine corner_check


end module remove_broken_corners