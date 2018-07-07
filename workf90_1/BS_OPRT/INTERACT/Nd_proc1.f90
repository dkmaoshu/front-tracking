module node_production_n_update_1
! This module produces new node cell for the merged critical cell
! and deals with the related update of discontinuity curves. The 
! node-cell information, such as the ordinary cell-average, the
! coordinates of the node, et al., and the discontinuity curve
! connecting the old and new node cells will be provided in the
! second version of the modules.

use variables_4_node_treat
! 'bs_var.f90'

implicit none

type(curve_plug), pointer :: tempp_old, tempp_posi, tempp_new
character*6 :: previous_merged
! Tell the merged critical cell of the previous plug insert. This
! information is used in collecting conserved quantity(ies).
character*6 :: side_deleted
! Tell the side from which the critical cell is going ot be deleted.

public  plugs_update
private	update_plug_rings, creat_new_plug, plug_update, &
        collect_conserved, tempp_old, tempp_new, tempp_posi, &
		previous_merged, side_deleted


contains


subroutine plugs_update
! This subroutine performces the production and update.

implicit none
! type(state) :: suum

nullify(tempp_new); nullify(tempp_old); nullify(tempp_posi)

! suum=sum
call creat_node(nd_nmb_n)

! call check_ring(1,'counter ')

call update_plug_rings

! call check_ring(1,'counter ')

end subroutine plugs_update


subroutine update_plug_rings
! This subroutine builds plug-ring for the new node cell and
! update the plug-ring for the old one.

implicit none
character*6 :: merged
logical :: if_participate
integer :: idx_c, idx_n
 
! type(state) :: suum 
! type(cv_plug_info) :: plugg
! logical :: xxx
! type(adss_plug) :: addss

! suum=sum

idx_c=tempp%address%idx; idx_n=tempp_n%address%idx

! Move the current curve-plug under concern to the new node cell 
! as the last plug in the ring.
previous_merged='zeroth'
tempp_posi=>ndd(nd_nmb_n)%cv_rg%eend
tempp_old=>tempp
side_deleted=tempp_old%plug%side_front
call plug_update(current_merged,'before','before')
head_mark=ndd(nd_nmb)%cv_rg%begin%next%address%idx
! plugg=tempp_old%plug
! call check_ring(2,'counter ')

! Move curve-plugs of the discontinuity curves behind the current
! discontinuity curve that also participate in the mergence.
if_participate=.false.
previous_merged=current_merged
call check_merge(if_participate,merged)
do while(if_participate)
 tempp_posi=>tempp_new
 side_deleted=tempp_old%plug%side_front
 call plug_update(merged,'before','before')
 if(.not.associated(ndd(nd_nmb)%cv_rg%begin%next)) exit
 head_mark=ndd(nd_nmb)%cv_rg%begin%next%address%idx
 call check_merge(if_participate,merged)
 if(tempp_old%address%idx.eq.idx_n) if_participate=.false.
 previous_merged=merged

!  call check_ring(2,'counter ')

end do
 tempp=>tempp_old; tempn=>tempp_new
! plugg=tempp_old%plug
! call check_ring(2,'counter ')

! call check_list_c(2,'up    ')

! plugg=tempp%plug

! Move the next to the current curve-plug to the new node cell as
! the first plug in the ring.
previous_merged=current_merged
tempp_posi=>ndd(nd_nmb_n)%cv_rg%begin%next
tempp_old=>tempp_n
side_deleted=tempp_old%plug%side_behind
call plug_update(next_merged,'before','after ')
tempp_n=>tempp_old; tempn_n=>tempp_new


! plugg=tempp_n%plug
! call check_ring(2,'counter ')

! Move curve-plugs of the discontinuity curves ahead of the next to
! the current discontinuity curve that also participate in the 
! mergence.
if(associated(ndd(nd_nmb)%cv_rg%begin%next)) then
 head_mark=ndd(nd_nmb)%cv_rg%begin%next%address%idx
 if_participate=.false.
 previous_merged=next_merged
 call check_merge(if_participate,merged)
 do while(if_participate)
  tempp_posi=>tempp_new
  side_deleted=tempp%plug%side_behind
  call plug_update(merged,'after ','after ')
  head_mark=ndd(nd_nmb)%cv_rg%begin%next%address%idx
  call check_merge(if_participate,merged)
  if(.not.associated(tempp_old)) exit
  previous_merged=merged
 end do
end if

tempp_n=>tempp_old; tempn_n=>tempp_new
! plugg=tempn_n%plug

if(.not.associated(ndd(nd_nmb)%cv_rg%begin%next)) then
 current_node_empty='yes'
end if

! call check_list_c(4,'down  ')
 
! plugg=tempp%plug

if(.not.associated(tempp_old)) nullify(tempp)

end subroutine update_plug_rings


subroutine plug_update(merged,insert_facing,delete_facing)
! This subroutine updates curve-plugs in node splitting.

implicit none
character*6, intent(in) :: merged
character*6, intent(in) :: insert_facing, delete_facing
! Tell whether the transfer of curve-plug goes in clock or 
! counter-clock direction.

type(cv_plug_info) :: plug, plug_new
character*6 :: side_insert

! type(state) :: suum
! type(cv_plug_info) :: plugg
! suum=ndd(2)%n_cell%or_state

plug=tempp_old%plug
select case(insert_facing)
 case('before')
  side_insert=plug%side_behind
 case('after ')
  side_insert=plug%side_front
 case default; print*, 'Something wrong.'; pause
end select

! Pick the plug-in information from the will-be updated plug.
call creat_new_plug(plug,merged,plug_new)
! Creat the corresponding new plug-in information for the new 
! node cell.
allocate(tempp_new)
tempp_new%plug=plug_new

! call check_ring(2,'counter ')

! Form the pointer pointing to the new plug. 
call insert(tempp_posi,tempp_new,insert_facing)
! Insert the new plug in the plug-ring of the new node cell.
 
! call check_ring(2,'counter ')

call collect_conserved(merged,side_deleted)
! Collect the conserved physical quantity(ies) involved in the
! update.
 
!  call check_list_c(1,'up    ')

! call check_ring(1,'clock   ')

! plugg=tempp_old%plug

call deletee(tempp_old,delete_facing)
! Delete the current curve_plug from the curve-ring of the
! old node cell.

! plugg=tempp_old%plug

end subroutine plug_update


subroutine creat_new_plug(plug_old,merged,plug_new)
! This subroutine creats a new curve plug from an old curve plug
! according to the obtained merge information.

implicit none
type(cv_plug_info), intent(in) :: plug_old
character*6, intent(in) :: merged
! Tell whether the head or next to the head critical participates
! the merge.
type(cv_plug_info), intent(out) :: plug_new

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
character*3 :: if_empty
type(state), dimension(-1:1,-1:1) :: u_f, uu_f, u_b, uu_b
integer :: i, j, x_dif, y_dif

call clean_up_plug(plug_new)
plug_new=plug_old

! The following determines the plug-in edge.
nullify(temp)
call locate_head(plug_old,merged,temp,if_empty)
select case(merged)
 case('zeroth')
 case('first ','second')
  g_cell=temp%g_cell; nullify(temp)
  if(plug_new%end_type.eq.'begin ') then
   plug_new%edge=g_cell%edge(2)
  else
   plug_new%edge=g_cell%edge(1)
  end if
 case default; print*, 'Something wrong!!!.'; pause
end select

! The following moves the local smooth in the old curve-plug to
! the new one.
do i=-1, 1; do j=-1, 1
 uu_f(i,j)=error_data; uu_b(i,j)=error_data 
 u_f(i,j)=plug_old%u_front(i,j); u_b(i,j)=plug_old%u_behind(i,j)
end do; end do
x_dif=x_nd_new-x_nd; y_dif=y_nd_new-y_nd
do i=-1, 1; do j=-1, 1
 if(iabs(i+x_dif).le.1.and.iabs(j+y_dif).le.1) then
  uu_f(i,j)=u_f(i+x_dif,j+y_dif)
  uu_b(i,j)=u_b(i+x_dif,j+y_dif)
 end if
end do; end do
do i=-1, 1; do j=-1, 1
 plug_new%u_front(i,j)=uu_f(i,j); plug_new%u_behind(i,j)=uu_b(i,j)
end do; end do

end subroutine creat_new_plug


subroutine collect_conserved(merged,side_deleted)
! This subroutine collect the conserved physical quantity(ies)
! involved in the update.

implicit none
character*6, intent(in) :: merged
! Tell whether the first or the second head critical cell is going
! to be merged.
character*6 , intent(in) :: side_deleted
! Tell whether the plug comes in from the fornt or the back. 

type(critical_cell), pointer :: temp, temp_next
type(cv_plug_info) :: plug
type(geo_info) :: g_cell, g_cell_next
type(phy_info) :: p_cell, p_cell_next
type(state) :: difference
character*6 :: posi
character*3 :: if_empty

nullify(temp); nullify(temp_next)
call clean_up_g_cell(g_cell); call clean_up_g_cell(g_cell_next)
call clean_up_p_cell(p_cell); call clean_up_p_cell(p_cell_next)
if_empty='no '; difference=error_data
plug=tempp_old%plug
call locate_head(plug,'first ',temp,if_empty)
g_cell=temp%g_cell; p_cell=temp%p_cell; if_empty='no '
call locate_head(plug,'second',temp_next,if_empty)
if(associated(temp_next)) then
 g_cell_next=temp_next%g_cell
 p_cell_next=temp_next%p_cell
end if
! Get geometrical and physical information of the head critical 
! cells of the dicontinuity curve whose plug will be updated.

select case(plug%end_type)
 case('begin '); posi='after '
 case('end   '); posi='before'
 case default;   call error_message
end select

select case(merged)
 case('zeroth')

 case('first ')
  sum=sum+p_cell%or_state
  if(previous_merged.eq.'first '.or.previous_merged.eq.'second') then
   select case(side_deleted)
    case('left  '); sum=sum-p_cell%l_state
    case('right '); sum=sum-p_cell%r_state
    case default;   call error_message
   end select 		 
  end if
  call remove(temp,posi,side_deleted,difference)

 case('second')
  call remove(temp,posi,side_deleted,difference)
  sum=sum+difference
  sum=sum+p_cell_next%or_state
  if(previous_merged.eq.'first '.or.previous_merged.eq.'second') then
   select case(side_deleted)
    case('left  '); sum=sum-p_cell_next%l_state
    case('right '); sum=sum-p_cell_next%r_state
    case default;   call error_message
   end select 		 
  end if
  call remove(temp,posi,side_deleted,difference)

end select

! The following updates the corresponding end of the curve.
select case(plug%end_type)
 case('begin '); cvv(plug%cv_nmb)%begin_end=nd_nmb_n
 case('end   '); cvv(plug%cv_nmb)%end_end=nd_nmb_n
 case default; call error_message
end select

! sum=sum-ndd(nd_nmb_n)%n_cell%or_state

end subroutine collect_conserved


subroutine check_merge(if_participate,merged)
! This subroutine checks a curve-plug to see if its head critical
! cell should participate the mergence.

implicit none
logical, intent(out) :: if_participate
character*6, intent(out) :: merged

type(critical_cell), pointer :: temp
type(cv_plug_info) :: plug
type(geo_info) :: g_cell
character*3 :: if_empty
integer :: i_c, j_c, i_n, j_n

if_participate=.false.; merged='      '
call clean_up_g_cell(g_cell)
plug=tempp_old%plug; nullify(temp)
call locate_head(plug,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
g_cell=temp%g_cell; nullify(temp)
if(if_empty.ne.'yes') then
 i_c=g_cell%x_idx; j_c=g_cell%y_idx
 if(i_c.eq.x_nd_new.and.j_c.eq.y_nd_new) then
  if_participate=.true.; merged='first '
 else
  if(plug%end_type.eq.'begin ') then
   call find_neighbor_cell(i_c,j_c,g_cell%edge(1),i_n,j_n)
  else
   call find_neighbor_cell(i_c,j_c,g_cell%edge(2),i_n,j_n)
  end if
  if(i_n.eq.x_nd_new.and.j_n.eq.y_nd_new) then
   if_participate=.true.; merged='zeroth'
  else
   call clean_up_g_cell(g_cell)
   call locate_head(plug,'second',temp,if_empty)
   if(associated(temp)) then
    g_cell=temp%g_cell; nullify(temp)
    i_c=g_cell%x_idx; j_c=g_cell%y_idx
    if(i_c.eq.x_nd_new.and.j_c.eq.y_nd_new) then
     if_participate=.true.; merged='second'
    end if
   end if
  end if
 end if
end if

end subroutine check_merge


end module node_production_n_update_1