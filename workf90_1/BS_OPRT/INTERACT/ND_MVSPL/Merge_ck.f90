module merge_checks
! This module contains two checks on whether there is merge of
! discontinuity curves.

use variables_4_node_treat
! bs_var.f90'

public  check_mergence
private check1, check2, check3, check4

! In the following code, when refering diagonal or neighboring
! cells the view is taken in the node cell under concern.


contains


subroutine check_mergence
! This subroutine implements the three checks of mergence.

implicit none

! call check_ring(1,'counter ')

if_merge='no '
call check1
! The first check on whether there is mergence of curves. 
if(if_merge.eq.'no ') call check2
! The second and third checks on whether there is mergence of 
! curves.
if(if_merge.eq.'no ') call check3

if(if_merge.eq.'no ') call check4
! The fourth check on whether there is mergence of curves. 

end subroutine check_mergence


subroutine check1
! This subroutine checks the case that the head critical cells of
! two neighboring discontinuity curves are in the same grid cell.

implicit none

type(critical_cell), pointer :: temp, temp_second
type(cv_plug_info) :: plug, plug_n
type(geo_info) :: g_cell, gn_cell, g_second
character*3 :: if_empty
integer :: i0, j0, i1, j1, num_c, num_n, diff_0, diff_1
real(8), dimension(2) :: dif_1, dif_2
logical :: corner_c, corner_n, if_right_direction
real(8) :: length

tempp_n=>tempp%next
plug=tempp%plug; plug_n=tempp_n%plug
if(cvv(plug%cv_nmb)%status.ne.'awake '.or. &
   cvv(plug_n%cv_nmb)%status.ne.'awake ') return

! Locate the pointers pointing to the head critical cells of the 
! curves under concern.
nullify(temp)
call locate_head(plug,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
g_cell=temp%g_cell; nullify(temp)
call locate_head(plug_n,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
gn_cell=temp%g_cell; nullify(temp)

! If the two head critical cell are not located in the same grid
! cell, check is done with a negative answer for 'if_merge'.
i0=g_cell%x_idx; j0=g_cell%y_idx
i1=gn_cell%x_idx; j1=gn_cell%y_idx
if(i0.ne.i1.or.j0.ne.j1) then
 if_merge='no '; return
end if

! Kick out the corner case.
corner_c=.false.; corner_n=.false.
nullify(temp_second)
diff_0=iabs(i0-x_nd)+iabs(j0-y_nd)
diff_1=iabs(i1-x_nd)+iabs(j1-y_nd)
if(diff_0.gt.1.and.g_cell%g_type.eq.'xy ') then
 call locate_head(plug,'second',temp_second,if_empty)
 if(associated(temp_second)) then
  g_second=temp_second%g_cell 
  if(iabs(g_second%x_idx-x_nd)+iabs(g_second%y_idx-y_nd).eq.1) then
   corner_c=.true.
  end if
 end if
end if
nullify(temp_second)
if(diff_1.gt.1.and.gn_cell%g_type.eq.'xy ') then
 call locate_head(plug_n,'second',temp_second,if_empty)
 if(associated(temp_second)) then
  g_second=temp_second%g_cell
  if(iabs(g_second%x_idx-x_nd)+iabs(g_second%y_idx-y_nd).eq.1) then
   corner_n=.true.
  end if
 end if
end if
if(corner_c.and.corner_n) return

! The case that both the head critical cell moved to a diagonal
! grid cell.

if(iabs(i0-x_nd)+iabs(j0-y_nd).gt.1) then
! The the following finds the x- and y-coordinates of the node
! point in the will-be node cell.
 
 call check_direction
! Only when the head critical cell moves counterclockwisely and the next head critical cell
! moves clockwisely into the diagonal grid cell the mergence will take place. This will be
! checked by the above subroutine 'check_direction'.

 if(if_right_direction)	then
  dif_1=0.d0; dif_2=0.d0
  call find_node(g_cell,gn_cell,dif_1,dif_2,xw_posi,yw_posi)
  if(.not.corner_c.and..not.corner_n) then
   if_merge='yes'; current_merged='first '; next_merged='first '
  else
   if(dabs(xw_posi).lt.0.5d0.and.dabs(yw_posi).lt.0.5d0) then
    if_merge='yes'; current_merged='first '; next_merged='first '
   end if
  end if
 else
  return 
 end if

end if
 
! The case that neither head critical cell moved. 

! The treatment of this case is just provisional; further developed 
! version of the algorithm will move this treatment to part that
! deals with interactions of discontinuity curves.

if(iabs(i0-x_nd)+iabs(j0-y_nd).eq.1) then
! Pick the numbers of the discontinuity positions of the two head
! critical cell that are on the edge of the node cell.
 if(neighbor_edge(g_cell%edge(1)).eq.plug%edge) then
  num_c=1
 else
  num_c=2
 end if
 if(neighbor_edge(gn_cell%edge(1)).eq.plug_n%edge) then
  num_n=1
 else
  num_n=2
 end if

! Check whether there is mergence.
 select case( plug%edge)
  case(1,2)
   if(gn_cell%dis_posi(num_n).lt.g_cell%dis_posi(num_c)+1.0d-8) then
    if_merge='yes'; current_merged='first '; next_merged='first '
   end if   
  case(3,4)
   if(gn_cell%dis_posi(num_n).gt.g_cell%dis_posi(num_c)) then
    if_merge='yes'; current_merged='first '; next_merged='first '
   end if   
 end select
end if

! If there is mergence, find the x- and y-coordinates of the node
! point in the will-be node cell.
dif_1=0.d0; dif_2=0.d0
if(if_merge.ne.'no ') then
 call compute_curve_length(g_cell,length)
 nullify(temp_second)
 if(length.lt.1.0d-2) then
  call locate_head(plug,'second',temp_second,if_empty)
  if(associated(temp_second)) then
   g_second=temp_second%g_cell
   dif_1(1)=dfloat(g_second%x_idx-g_cell%x_idx)
   dif_1(2)=dfloat(g_second%y_idx-g_cell%y_idx)
   g_cell=g_second
  else
   call error_message
  end if
 end if
 call compute_curve_length(gn_cell,length)
 nullify(temp_second)
 if(length.lt.1.0d-2) then
  call locate_head(plug_n,'second',temp_second,if_empty)
  if(associated(temp_second)) then
   g_second=temp_second%g_cell
   dif_2(1)=dfloat(g_second%x_idx-gn_cell%x_idx)
   dif_2(2)=dfloat(g_second%y_idx-gn_cell%y_idx)
   gn_cell=g_second
  else
   call error_message
  end if
 end if
 call find_node(g_cell,gn_cell,dif_1,dif_2,xw_posi,yw_posi)
 if(dabs(xw_posi).gt.0.5d0.or.dabs(yw_posi).gt.0.5d0) then
  if(dabs(xw_posi).gt.0.5d0) xw_posi=dsign(0.4999999d0,xw_posi)
  if(dabs(yw_posi).gt.0.5d0) yw_posi=dsign(0.4999999d0,yw_posi)
 ! xw_posi=0.d0; yw_posi=0.d0
 end if
 x_nd_new=i0; y_nd_new=j0
end if


contains


subroutine check_direction

implicit none

character*3 :: position
! Indicating the position of the diagonal grid cell.

if_right_direction=.false.

if(i0.eq.x_nd) call error_message
if(j0.eq.y_nd) call error_message
if(i0.lt.x_nd) then
 if(j0.lt.y_nd) then
  position='s-w'
 else
  position='n-w'
 end if
else
 if(j0.lt.y_nd) then
  position='s-e'
 else
  position='n-e'
 end if
end if

select case(position)
 case('n-e')
  select case(plug%edge)
   case(2); if_right_direction=.true.
   case(1,4); call error_message
  end select
 case('n-w')
  select case(plug%edge)
   case(3); if_right_direction=.true.
   case(1,2); call error_message
  end select
 case('s-e')
  select case(plug%edge)
   case(1); if_right_direction=.true.
   case(3,4); call error_message
  end select
 case('s-w')
  select case(plug%edge)
   case(4); if_right_direction=.true.
   case(2,3); call error_message
  end select
end select
  
end subroutine check_direction


end subroutine check1


subroutine check2
! This subroutine checks the case that one of the head critical 
! cells moved, which is now in a diagonal grid cell, and the 
! othetr stayed, which is now in a neighboring grid cell. The 
! crossing of the two curves could happen either in the staying
! head critical cell or in the critical cell next to the staying
! one on the curve.

implicit none

type(critical_cell), pointer :: temp, temp_second
type(cv_plug_info) :: plug, plug_n
type(geo_info) :: g_cell, gn_cell, gm_cell, gs_cell, g_second
character*3 :: if_empty
integer :: i0, j0, i1, j1, dif0, dif1, x_staying, y_staying
real(8), dimension(2) :: dif_1, dif_2
character*8 :: staying, corner_1, corner_2, corner_m
real(8) :: length

tempp_n=>tempp%next
plug=tempp%plug; plug_n=tempp_n%plug
if(cvv(plug%cv_nmb)%status.ne.'awake '.or. &
   cvv(plug_n%cv_nmb)%status.ne.'awake ') return

! Locate the head critical cells of the curves under concern.
nullify(temp)
call locate_head(plug,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
g_cell=temp%g_cell; nullify(temp)
call locate_head(plug_n,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
gn_cell=temp%g_cell; nullify(temp)

i0=g_cell%x_idx; j0=g_cell%y_idx
i1=gn_cell%x_idx; j1=gn_cell%y_idx
dif0=iabs(i0-x_nd)+iabs(j0-y_nd)
dif1=iabs(i1-x_nd)+iabs(j1-y_nd)
if(dif0.gt.1.and.dif1.gt.1) return
if(dif0.le.1.and.dif1.le.1) return
if(iabs(i1-i0).gt.1.or.iabs(j1-j0).gt.1) return

! Kick out the corner case.
corner_1='no '; corner_2='no '
nullify(temp_second)
if(dif0.gt.1.and.g_cell%g_type.eq.'xy ') then
 call locate_head(plug,'second',temp_second,if_empty)
 if(associated(temp_second)) then
  g_second=temp_second%g_cell
  if(iabs(g_second%x_idx-x_nd)+iabs(g_second%y_idx-y_nd).eq.1) corner_1='yes'
 end if
end if
nullify(temp_second)
if(dif1.gt.1.and.gn_cell%g_type.eq.'xy ') then
 call locate_head(plug_n,'second',temp_second,if_empty)
 if(associated(temp_second)) then
  g_second=temp_second%g_cell
  if(iabs(g_second%x_idx-x_nd)+iabs(g_second%y_idx-y_nd).eq.1) corner_2='yes'
 endif
end if

! Identify the head critical cell that moved and the one that stays.
if(dif0.gt.1) then
 gm_cell=g_cell; gs_cell=gn_cell; staying='next    '; corner_m=corner_1
else
 if(dif1.gt.1) then
  gm_cell=gn_cell; gs_cell=g_cell; staying='current '; corner_m=corner_2
 else
  print*, 'There must be something wrong'; pause
 end if
end if

x_staying=gs_cell%x_idx; y_staying=gs_cell%y_idx
call compute_curve_length(gs_cell,length)
if(length.lt.1.0d-2) then
 nullify(temp_second)
 select case(staying)
  case('current ') 
   call locate_head(plug,'second',temp_second,if_empty)
   if(associated(temp_second)) then
    g_second=temp_second%g_cell
   else
    return
   end if
  case('next    ')
   call locate_head(plug_n,'second', temp_second,if_empty)
   if(associated(temp_second)) then
    g_second=temp_second%g_cell
   else
    return
   end if	 		   	 
 end select
 gs_cell=g_second; nullify(temp_second); call clean_up_g_cell(g_second)
end if
dif_1(1)=dfloat(gs_cell%x_idx-x_staying)
dif_1(2)=dfloat(gs_cell%y_idx-y_staying)


call compute_curve_length(gm_cell,length)
if(length.lt.1.0d-2) then
 nullify(temp_second)
 select case(staying)
  case('current ') 
   call locate_head(plug_n,'second',temp_second,if_empty)
   if(associated(temp_second)) then
    g_second=temp_second%g_cell
   else
    return
   end if
  case('next    ')
   call locate_head(plug,'second', temp_second,if_empty)
   if(associated(temp_second)) then
    g_second=temp_second%g_cell
   else
    return
   end if	 		   	 
 end select
 gm_cell=g_second; nullify(temp_second); call clean_up_g_cell(g_second)
end if
dif_2(1)=dfloat(gm_cell%x_idx-x_staying)
dif_2(2)=dfloat(gm_cell%y_idx-y_staying)

! Find the intersection point of the curve segment in the staying
! critical cell and the extension of the curve segment in the
! moving critical cell and then determine the merge situation.
if(corner_m.eq.'no ') then
 call find_node(gs_cell,gm_cell,dif_1,dif_2,xw_posi,yw_posi)
 if(dabs(xw_posi).lt.0.5d0.and.dabs(yw_posi).lt.0.5d0) then
  if_merge='yes'
  if(staying.eq.'current ') then
   current_merged='first '; next_merged='zeroth'
  else
   current_merged='zeroth'; next_merged='first '
  end if
 end if
end if

! Find the crossing point of the curve segment in the critical 
! cell next to the staying one and the extension of the curve 
! segment in the moving critical cell and then determine the
! merge situation.

if(if_merge.eq.'no ') then
 if(staying.eq.'current ') then
  call locate_head(plug,'second',temp,if_empty)
 else
  call locate_head(plug_n,'second',temp,if_empty)
 end if
 if(.not.associated(temp)) return
 gs_cell=temp%g_cell
 if(gm_cell%x_idx.ne.gs_cell%x_idx.or.gm_cell%y_idx.ne.gs_cell &
	%y_idx) return
 dif_1=0.d0; dif_2=0.d0
 call find_node(gs_cell,gm_cell,dif_1,dif_2,xw_posi,yw_posi)
 if(dabs(xw_posi).lt.0.5d0.and.dabs(yw_posi).lt.0.5d0) then
  if_merge='yes'
  if(staying.eq.'current ') then
   current_merged='second'; next_merged='first '
  else
   current_merged='first '; next_merged='second'
  end if
 end if
end if

if(if_merge.eq.'yes') then
 x_nd_new=gs_cell%x_idx; y_nd_new=gs_cell%y_idx
end if

end subroutine check2


subroutine check3

implicit none

type(critical_cell), pointer :: temp, temp_second
type(cv_plug_info) :: plug, plug_n
type(geo_info) :: g_cell, gn_cell, gm_cell, gs_cell,  &
                  g_second, gn_second
character*3 :: if_empty
integer :: i0, j0, i1, j1
real(8), dimension(2) :: dif_1, dif_2
character*8 :: staying
logical :: corner_1, corner_2

tempp_n=>tempp%next
plug=tempp%plug; plug_n=tempp_n%plug
if(cvv(plug%cv_nmb)%status.ne.'awake '.or. &
   cvv(plug_n%cv_nmb)%status.ne.'awake ') return

! Locate the head critical cells of the curves under concern.
nullify(temp)
call locate_head(plug,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
g_cell=temp%g_cell; nullify(temp)
call locate_head(plug_n,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
gn_cell=temp%g_cell; nullify(temp)

i0=g_cell%x_idx; j0=g_cell%y_idx
i1=gn_cell%x_idx; j1=gn_cell%y_idx
if(i0.ne.i1.or.j0.ne.j1) return
if(iabs(i0-x_nd).ne.1.or.iabs(j0-y_nd).ne.1) return

! Kick out the corner case.
corner_1=.false.; corner_2=.false.
call clean_up_g_cell(g_second); call clean_up_g_cell(gn_second)
nullify(temp_second)
if(g_cell%g_type.eq.'xy ') then
 call locate_head(plug,'second',temp_second,if_empty)
 if(associated(temp_second)) then
  g_second=temp_second%g_cell
  if(iabs(g_second%x_idx-x_nd)+iabs(g_second%y_idx-y_nd).eq.1) then
   corner_1=.true.
  end if
 end if
end if
nullify(temp_second)
if(gn_cell%g_type.eq.'xy ') then
 call locate_head(plug_n,'second',temp_second,if_empty)
 if(associated(temp_second)) then
  gn_second=temp_second%g_cell
  if(iabs(gn_second%x_idx-x_nd)+iabs(gn_second%y_idx-y_nd).eq.1) then
   corner_2=.true.
  end if
 endif
end if
if(corner_1.and.corner_2) return
if(.not.corner_1.and..not.corner_2) call error_message

! Identify the head critical cell that moved and the one that stays.
if(.not.corner_1) then
 gm_cell=g_cell; gs_cell=gn_second; staying='next    '
else
 gm_cell=gn_cell; gs_cell=g_second; staying='current '
end if

! Find the intersection point of the curve segment in the staying
! critical cell and the extension of the curve segment in the
! moving critical cell and then determine the merge situation.
dif_1=0.d0
dif_2(1)=gm_cell%x_idx-gs_cell%x_idx
dif_2(2)=gm_cell%y_idx-gs_cell%y_idx
call find_node(gs_cell,gm_cell,dif_1,dif_2,xw_posi,yw_posi)
if(dabs(xw_posi).lt.0.5d0.and.dabs(yw_posi).lt.0.5d0) then
 x_nd_new=gs_cell%x_idx; y_nd_new=gs_cell%y_idx
 if_merge='yes'
 if(staying.eq.'current ') then
  current_merged='second'; next_merged='zeroth'
 else
  current_merged='zeroth'; next_merged='second'
 end if
end if

end subroutine check3


subroutine check4
! This subroutine checks the case that both the head critical  
! cells moved.

implicit none

type(critical_cell), pointer :: temp, temp_second
type(cv_plug_info) :: plug, plug_n
type(geo_info) :: g_cell, gn_cell, g_second
character*3 :: if_empty
integer :: i0, j0, i1, j1, edge_c0, edge_c1, dif0, dif1
integer :: ii0, jj0, ii1, jj1
real(8), dimension(2) :: dif_1,dif_2
real(8) :: posi1, posi2
character*10 :: error_msg

tempp_n=>tempp%next
plug=tempp%plug; plug_n=tempp_n%plug
if(cvv(plug%cv_nmb)%status.ne.'awake '.or. &
   cvv(plug_n%cv_nmb)%status.ne.'awake ') return

! Locate the head critical cells of the curves under concern.
nullify(temp)
call locate_head(plug,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
g_cell=temp%g_cell; nullify(temp)
call locate_head(plug_n,'first ',temp,if_empty)
if(if_empty.eq.'yes') return
gn_cell=temp%g_cell; nullify(temp)

i0=g_cell%x_idx; j0=g_cell%y_idx
i1=gn_cell%x_idx; j1=gn_cell%y_idx
dif0=iabs(i0-x_nd)+iabs(j0-y_nd)
dif1=iabs(i1-x_nd)+iabs(j1-y_nd)
if(dif0.le.1.or.dif1.le.1) return

! Kick out the corner case.
nullify(temp_second)
if(dif0.gt.1.and.g_cell%g_type.eq.'xy ') then
 call locate_head(plug,'second',temp_second,if_empty)
 if(associated(temp_second)) then
  g_second=temp_second%g_cell
  if(iabs(g_second%x_idx-x_nd)+iabs(g_second%y_idx-y_nd).eq.1) return
 end if
end if
nullify(temp_second)
if(dif1.gt.1.and.gn_cell%g_type.eq.'xy ') then
 call locate_head(plug_n,'second',temp_second,if_empty)
 if(associated(temp_second)) then
  g_second=temp_second%g_cell
  if(iabs(g_second%x_idx-x_nd)+iabs(g_second%y_idx-y_nd).eq.1) return
 end if
end if

! Find the cell edge of the node cell adjacent to the grid cell in
! between for each of the head critical cells under concern.
call find_edgec(g_cell,plug,ii0,jj0,edge_c0)
call find_edgec(gn_cell,plug_n,ii1,jj1,edge_c1)
if(edge_c0.ne.edge_c1.or.ii0.ne.ii1.or.jj0.ne.jj1) return

! Find the intersection point of the extension of the curve
! segment in each head critical cell and the common edge.
call find_plug_posi(g_cell,edge_c0,posi1,error_msg)
if(error_msg.eq.'parallel  '.or.dabs(posi1).ge.0.5d0) return
call find_plug_posi(gn_cell,edge_c1,posi2,error_msg)
if(error_msg.eq.'parallel  '.or.dabs(posi2).ge.0.5d0) return
select case(edge_c1)
 case(1,2) 
  if(posi1.gt.posi2) then
   if_merge='yes'; current_merged='zeroth'; next_merged='zeroth'
  end if
 case(3,4)
  if(posi1.lt.posi2) then
   if_merge='yes'; current_merged='zeroth'; next_merged='zeroth'
  end if
end select

if(if_merge.ne.'no ') then   
 dif_1(1)=g_cell%x_idx-ii0; dif_1(2)=g_cell%y_idx-jj0
 dif_2(1)=gn_cell%x_idx-ii1; dif_2(2)=gn_cell%y_idx-jj1
 call find_node(g_cell,gn_cell,dif_1,dif_2,xw_posi,yw_posi)
 if(dabs(xw_posi).gt.0.5d0.or.dabs(yw_posi).gt.0.5d0) then
  xw_posi=0.5d0; yw_posi=0.5d0
 end if
end if

if(if_merge.eq.'yes') then
 x_nd_new=ii0; y_nd_new=jj0
end if

end subroutine check4


end module merge_checks
