module fix_broken_connection
! This module fixes broken connections between the old nodel cell
! and discontinuity curves supposed to plug-in it.

use tools_4_connecting_n_triple
! 'tools_ct.f90'

implicit none

type(cv_plug_info) :: plug, pl_new
type(adss_plug) :: adss
type(critical_cell), pointer :: temp, t_new, temp_second
type(geo_info) :: g_cell, g_cell_s, g_new
type(phy_info) :: p_cell_s, p_new
character*3 :: if_empty
real(8) :: area
type(state) :: dif_cov
! The difference in conservation caused by the change of physical
! state in the grid cell.
integer :: single, st_edge,ed_edge, edd_edge, x_idx_s, y_idx_s
real(8), dimension(2) :: st_posi, ed_posi

public connecting
private	new_geo_cell_n_area, new_geo_cell_4_yawn, awake_curve, &
        yawn_curve, plug, pl_new, adss, temp, t_new, g_cell, &
		g_new, p_new, if_empty, area, dif_cov, single, st_edge, &
		ed_edge, st_posi, ed_posi, physical_states, yawn_update, &
        update_near_node, g_cell_s, temp_second, x_idx_s, y_idx_s, p_cell_s


contains


subroutine connecting

! This subroutine performs the connecting job

implicit none

integer :: i0, j0, i1, j1, end_edge
type(node_info) :: n_cell

! type(cv_plug_info) :: plugg
! type(adss_plug) :: adsss

! call check_list_c(5,'up    ')

tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)

 plug=tempp%plug; adss=tempp%address; call clean_up_plug(pl_new)
 call clean_up_g_cell(g_cell); call clean_up_p_cell(p_new)
 nullify(temp)
! Initialization.
 
 call locate_head(plug,'first ',temp,if_empty)
 if(associated(temp)) then
  g_cell=temp%g_cell
  if(iabs(x_nd-g_cell%x_idx)+iabs(y_nd-g_cell%y_idx).gt.1) then
! If broken of connection is found, connecting job sould be
! carried out.
   call awake_curve
  else
   select case(plug%end_type)
    case('begin '); end_edge=neighbor_edge(g_cell%edge(1))
    case('end   '); end_edge=neighbor_edge(g_cell%edge(2))
	case default; call error_message
   end select
   if(end_edge.ne.plug%edge) call error_message
! There must be something wrong if a disjoint is found.
  end if
 else
! For an empty 'yawn' discontinuity curve the whether-broken is
! checked on observation of locations of the two node cells the 
! discontinuity curve plugs in. 
  n_cell=ndd(cvv(plug%cv_nmb)%begin_end)%n_cell
  i0=n_cell%x_idx; j0=n_cell%y_idx
  n_cell=ndd(cvv(plug%cv_nmb)%end_end)%n_cell
  i1=n_cell%x_idx; j1=n_cell%y_idx

!   plugg=tempp%plug; adsss=tempp%address

  if(iabs(i1-i0)+iabs(j1-j0).gt.1) call yawn_curve
 end if

!  call check_ring(1,'clock   ')
!  call check_list_c(6,'up    ')

 if(associated(tempp)) then
  tempp=>tempp%next  
  head_switch=(tempp%address%idx.ne.head_mark)
 else
  head_switch=.false.
 end if
end do

end subroutine connecting


subroutine awake_curve
! This subroutine fixes connections for curves in 'awake' status 
! with the node cell.

implicit none

! logical :: xxx

! The following determines the indexes of the new critical cell.
g_new%x_idx=g_cell%x_idx; g_new%y_idx=g_cell%y_idx
select case(plug%edge)
 case(1,3); g_new%y_idx=y_nd
 case(2,4); g_new%x_idx=x_nd
 case default; print*, 'Something wrong!!!'; pause
end select
  
! The following forms the new head critical cell to make-up the
! connection.
! Compute geometrical information of the new critical cell first.
! Set the start and end information.
ed_posi(1)=ndd(nd_nmb)%n_cell%x_posi+x_nd-g_new%x_idx
ed_posi(2)=ndd(nd_nmb)%n_cell%y_posi+y_nd-g_new%y_idx
select case(plug%end_type)
 case('begin ')
  call pick_point(g_cell,1,st_posi); st_edge=g_cell%edge(1)
 case('end   ')
  call pick_point(g_cell,2,st_posi); st_edge=g_cell%edge(2)
 case default; print*, 'Something is wrong!!!'; pause
end select
st_edge=neighbor_edge(st_edge)
st_posi(1)=st_posi(1)+g_cell%x_idx-g_new%x_idx
st_posi(2)=st_posi(2)+g_cell%y_idx-g_new%y_idx
select case(st_edge)
 case(1,3)
  if(g_cell%x_idx.gt.x_nd) ed_edge=4
  if(g_cell%x_idx.lt.x_nd) ed_edge=2
 case(2,4)
  if(g_cell%y_idx.gt.y_nd) ed_edge=1
  if(g_cell%y_idx.lt.y_nd) ed_edge=3
end select

call new_geo_cell_n_area(st_posi,st_edge,ed_posi,ed_edge, &
                         plug%end_type,g_new,area,edd_edge)
if(g_new%g_type.ne.'xy ') then
 print*, 'Something is wrong!!!'; pause
end if

! Compute physical information of the new critical cell then.
call physical_states(plug,g_new,area,p_new)
p_new%wv_nb=temp%p_cell%wv_nb

pl_new=plug
select case(pl_new%end_type)
 case('begin '); pl_new%edge=neighbor_edge(g_new%edge(1))
 case('end   '); pl_new%edge=neighbor_edge(g_new%edge(2))
end select

! The following updates the numerical solution.
! First, update the curve-plug
tempp%plug=pl_new
! Next plug in the new critical cell to complete the connection.
allocate(t_new); t_new%g_cell=g_new; t_new%p_cell=p_new
select case(pl_new%end_type)
 case('begin '); call insert(temp,t_new,'before')
 case('end   '); call insert(temp,t_new,'after ')
end select

! Finally, update the solution in smooth, update the grid-map and
! preserve conservation.
call update_near_node(t_new,tempp)
call find_single(g_new,single)
select case(g_new%point(single))
 case('left  '); dif_cov=p_new%or_state-p_new%r_state
 case('right '); dif_cov=p_new%or_state-p_new%l_state
 case default; print*, 'Something is wrong!!!'; pause
end select
ndd(nd_nmb)%n_cell%or_state=ndd(nd_nmb)%n_cell%or_state-dif_cov

end subroutine awake_curve


subroutine yawn_curve
! This subroutine fixes connections for curves in 'yawn' status
! with the node cell.

implicit none

type(node_info) :: n_b, n_e
real(8), dimension(2) :: posi
type(state) :: l_state, r_state
type(phy_info), dimension(wave_number) :: p_new
type(phy_info), dimension(wave_number,1) :: pp_new
type(geo_info), dimension(1) :: gg_new
character*8 :: dir, insert_order
integer :: i, dis_number
character*12, dimension(wave_number) :: wave_type
logical :: if_first

! type(cv_plug_info) :: plugg
! type(adss_plug) :: adsss

! plugg=tempp%plug; adsss=tempp%address

if(cvv(plug%cv_nmb)%begin_end.eq.nd_nmb) then
 n_b=ndd(cvv(plug%cv_nmb)%end_end)%n_cell
 n_e=ndd(cvv(plug%cv_nmb)%begin_end)%n_cell
end if
if(cvv(plug%cv_nmb)%end_end.eq.nd_nmb) then
 n_b=ndd(cvv(plug%cv_nmb)%begin_end)%n_cell
 n_e=ndd(cvv(plug%cv_nmb)%end_end)%n_cell
end if

! The following finds out geometrical information for the new 
! critical cell.

! Set indices for the new critical cell.
select case(plug%edge)
 case(1); g_new%x_idx=n_e%x_idx; g_new%y_idx=n_e%y_idx-1
 case(2); g_new%x_idx=n_e%x_idx+1; g_new%y_idx=n_e%y_idx
 case(3); g_new%x_idx=n_e%x_idx; g_new%y_idx=n_e%y_idx+1
 case(4); g_new%x_idx=n_e%x_idx-1; g_new%y_idx=n_e%y_idx
 case default; print*, 'Something is wrong!!!'; pause
end select
! Set start- and end-edges
st_edge=neighbor_edge(plug%edge)
if(n_b%x_idx.gt.g_new%x_idx) ed_edge=2
if(n_b%x_idx.lt.g_new%x_idx) ed_edge=4
if(n_b%y_idx.gt.g_new%y_idx) ed_edge=1
if(n_b%y_idx.lt.g_new%y_idx) ed_edge=3
! Set possible positions
posi=0.4999999

call new_geo_cell_4_yawn(st_edge,ed_edge,plug%end_type,posi,g_new)

! The following sets the physical information for the new 
! critical cell.

! Set left and right states and the ring direction of the
! corresponding plugs.
select case(plug%end_type)
 case('begin ');
  l_state=plug%u_front(g_new%x_idx-x_nd,g_new%y_idx-y_nd)
  r_state=plug%u_behind(g_new%x_idx-x_nd,g_new%y_idx-y_nd)
  dir='clock   '
 case('end   ')
  l_state=plug%u_behind(g_new%x_idx-x_nd,g_new%y_idx-y_nd)
  r_state=plug%u_front(g_new%x_idx-x_nd,g_new%y_idx-y_nd)
  dir='counter '
end select
if_first=.true.

call physical_separation(l_state,r_state,g_new,dir,p_new, &
     wave_type,if_first,dis_number)

! p_new%l_state=l_state; p_new%r_state=r_state

do i=1, wave_number; pp_new(i,1)=p_new(i); end do
gg_new(1)=g_new

! The following updates the numerical solution.
 
! plugg=tempp%plug; adsss=tempp%address

call yawn_update(gg_new,pp_new,1,insert_order)

end subroutine yawn_curve



end module fix_broken_connection