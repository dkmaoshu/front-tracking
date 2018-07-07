module fix_triple_nodes
! This module provides the treatment (or fix) to triple node 
! either with 'awake'-status or 'yawn'-status outgoing 
! discontinuity curves.

use tools_4_connecting_n_triple
! 'tools_ct.f90'

use tools_used_in_triple
! 'trip_bg.f90'

implicit none

type cr_pointer
 type(critical_cell), pointer :: temp
 type(curve_plug), pointer :: temppx
end type cr_pointer

type(cr_pointer), dimension(:), allocatable :: crpts
type(curve_plug), pointer :: tempp_outgoing
type(cv_plug_info) :: plug, pl_new
type(adss_plug) :: adss
type(critical_cell), pointer :: temp, t_new
type(geo_info) :: g_cell, g_new
type(phy_info) :: p_new
character*3 :: if_empty
real(8) :: area
type(state) :: dif_cov
! The difference in conservation caused by the change of physical
! state in the grid cell.
type(state), dimension(:,:), allocatable :: dist_cov
! The difference in conservation distributed among critical cells
! on 'outgoing' discontinuity curves.
integer :: single, st_edge, ed_edge, cv_outgoing, edd_edge
real(8), dimension(2) :: st_posi, ed_posi
type(state), dimension(:), allocatable :: dist_conv
integer :: number_of_outgoing
! The number of outgoing discontinuity curves.

type(cr_pointer) :: crpts_i

public awake_curves, yawn_curve, tempp, tempp_outgoing
private plug, pl_new, adss,temp,t_new, g_cell, g_new, p_new, &
        if_empty, area, dif_cov, single, st_edge, ed_edge, &
		st_posi, ed_posi,new_geo_cell_n_area, physical_states, &
        new_geo_cell_4_yawn, yawn_update, cv_outgoing, edd_edge, &
		insert_n_update, set_start, decompose, distribution, &
		dist_cov, dist_conv, crpts


contains


subroutine awake_curves
! This subroutine deals with the case of 'awake'-status outgoing
! discontinuity curves.

implicit none

!type array_plug_pointer
! type(curve_plug), pointer :: tempp
!end type

type(curve_plug), pointer :: temppx, temppx_b, tpp_new, tempp_swop
!type(array_plug_pointer), dimension(:), allocatable :: temppx
type(cv_plug_info) :: pplug1, pplug2
integer :: i, j, i0, j0, i1, j1, i_dif, j_dif, nd_mem
character*3 :: dis_old, old_new
type(state) :: or_state
logical :: need_end_type_change
character*6 :: original_end_type, new_end_type

type(state), dimension(:), allocatable :: xxx

! type(adss_plug) :: address

! type(cv_plug_info) :: pppl

nd_mem=nd_nmb
pplug1=tempp_outgoing%plug
allocate(crpts(wave_number))
!allocate(temppx(wave_number+1))
cv_outgoing=tempp_outgoing%plug%cv_nmb
i0=ndd(nd_nmb)%n_cell%x_idx; j0=ndd(nd_nmb)%n_cell%y_idx
do i=1, wave_number
 nullify(crpts(i)%temp)
end do
call find_other_end(tempp_outgoing,temppx_b)
tempp=>temppx_b

number_of_outgoing=0
temppx=>tempp_outgoing%next
pplug1=temppx%plug; pplug2=temppx_b%plug
do while(temppx%plug%inout.eq.'outgoing')
 number_of_outgoing=number_of_outgoing+1
 plug=temppx%plug

! Plug the discontinuity curve into the new node cell.
 allocate(tpp_new)
 tpp_new%address%nd_nmb=nd_nmb_n
 pl_new=plug
 i_dif=x_nd_new-x_nd; j_dif=y_nd_new-y_nd
 do i=-1,1; do j=-1,1
  pl_new%u_front(i,j)=plug%u_front(i+i_dif,j+j_dif)
  pl_new%u_behind(i,j)=plug%u_behind(i+i_dif,j+j_dif)
 end do; end do
 pl_new%inout='        '
 tpp_new%plug=pl_new
 crpts(number_of_outgoing)%temppx=>tpp_new
 call insert(temppx_b,tpp_new,'after ')

! To see if we need to change the end_type of the outgoing plug.
 need_end_type_change=.false. 
 if(tempp%plug%end_type.ne.temppx%plug%end_type) then
  need_end_type_change=.true.
  original_end_type=tempp%plug%end_type
  new_end_type=temppx%plug%end_type
 end if
!  pppl=crpts(number_of_outgoing)%temppx%plug

 dis_old='no '; old_new='no '

! Inspect whether there should be a critical cell between the
! old node cell and the head of the discontinuity curve.
 i1=error_index; j1=error_index
 call locate_head(plug,'first ',temp,if_empty)
 if(if_empty.eq.'yes') then
  print*, 'Something is wrong!!!'; pause
 end if
 g_cell=temp%g_cell
 i1=g_cell%x_idx; j1=g_cell%y_idx
 if(iabs(i1-i0)+iabs(j1-j0).gt.1) dis_old='yes'
! Inspect whether there should be a critical cell between the
! old node cell and the new node cell.
 i1=ndd(nd_nmb_n)%n_cell%x_idx; j1=ndd(nd_nmb_n)%n_cell%y_idx
 if(iabs(i1-i0)+iabs(j1-j0).gt.1) old_new='yes'
 
! If 'dis_old' takes 'yes', produce a critical cell to extend
! the discontinuity in the cell.
 if(dis_old.eq.'yes') then  
! The following determines the indexes of the new critical cell.
  g_new%x_idx=g_cell%x_idx; g_new%y_idx=g_cell%y_idx
  select case(plug%edge)
   case(1,3); g_new%y_idx=y_nd
   case(2,4); g_new%x_idx=x_nd
   case default; print*, 'Something wrong!!!'; pause
  end select

! The following forms the new critical cell.
! Compute geometrical information of the new critical cell first.
! Set the start and end information.
  call set_start
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
! The following insert the critical cell into the curve and 
! update the numerical solution.
  if(need_end_type_change) tempp%plug%end_type=new_end_type

  tempp_swop=>tempp
  tempp=>tpp_new
  call insert_n_update('no ')
  tempp=>tempp_swop

  if(need_end_type_change) tempp%plug%end_type=original_end_type
   
  g_cell=g_new; call clean_up_g_cell(g_new)
 end if

!  pppl=crpts(1)%temppx%plug

! Form the critical cell in old node cell  
! The following determines the indexes of the new critical cell.
 g_new%y_idx=y_nd; g_new%x_idx=x_nd
! The following forms the new critical cell.
! Compute geometrical information of the new critical cell first.
! Set the start and end information.
 call set_start
 if(dis_old.eq.'yes'.and.edd_edge.gt.0) ed_edge=edd_edge
 call new_geo_cell_n_area(st_posi,st_edge,ed_posi,ed_edge, &
                          plug%end_type,g_new,area,edd_edge)
 select case(plug%end_type)
  case('begin '); plug%edge=neighbor_edge(g_new%edge(1))
  case('end   '); plug%edge=neighbor_edge(g_new%edge(2))
 end select
! Compute physical information of the new critical cell then.
 call physical_states(plug,g_new,area,p_new)
 p_new%wv_nb=temp%p_cell%wv_nb
! The following updates the numerical solution.
! Plug in the new critical cell to complete the connection.
 if(need_end_type_change) tempp%plug%end_type=new_end_type

 tempp_swop=>tempp
 tempp=>tpp_new
 call insert_n_update('yes')
 tempp=>tempp_swop

!  call check_list_c(5,'down  ')

 if(need_end_type_change) tempp%plug%end_type=original_end_type
 g_cell=g_new; call clean_up_g_cell(g_new)
! This critical cell should be pointed so that it can be picked 
! up later in preserving conservation.
 crpts(number_of_outgoing)%temp=>t_new

!  pppl=crpts(1)%temppx%plug
!  call check_list_c(3,'up    ')
 
! If 'old_new' takes 'yes', produce a critical cell to extend
! the discontinuity in the cell.
 if(old_new.eq.'yes') then
  g_new%x_idx=g_cell%x_idx; g_new%y_idx=g_cell%y_idx
  select case(plug%edge)
   case(1,3); g_new%y_idx=y_nd_new
   case(2,4); g_new%x_idx=x_nd_new
   case default; print*, 'Something wrong!!!'; pause
  end select

! The following forms the new critical cell.
! Compute geometrical information of the new critical cell first.
! Set the start and end information.
  call set_start
  select case(st_edge)
   case(1,3)
    if(g_cell%x_idx.gt.x_nd_new) ed_edge=4
    if(g_cell%x_idx.lt.x_nd_new) ed_edge=2
   case(2,4)
    if(g_cell%y_idx.gt.y_nd_new) ed_edge=1
    if(g_cell%y_idx.lt.y_nd_new) ed_edge=3
  end select
  call new_geo_cell_n_area(st_posi,st_edge,ed_posi,ed_edge, &
                           plug%end_type,g_new,area,edd_edge)
  select case(plug%end_type)
   case('begin '); plug%edge=neighbor_edge(g_new%edge(1))
   case('end   '); plug%edge=neighbor_edge(g_new%edge(2))
  end select
  if(g_new%g_type.ne.'xy ') then
   print*, 'Something is wrong!!!'; pause
  end if
! Compute physical information of the new critical cell then.
  call physical_states(plug,g_new,area,p_new)
  p_new%wv_nb=temp%p_cell%wv_nb
! The following insert the critical cell into the curve and 
! update the numerical solution.

!   call check_list_c(4,'up    ')
  if(need_end_type_change) tempp%plug%end_type=new_end_type

  tempp_swop=>tempp
  tempp=>tpp_new
  call insert_n_update('no ')
  tempp=>tempp_swop

  if(need_end_type_change) tempp%plug%end_type=original_end_type

!   pppl=temppx_b%plug
!   pppl=tempp%plug

!   call check_list_c(5,'down  ')

  g_cell=g_new; call clean_up_g_cell(g_new)
 end if

! Change the new plug cell-edge.
 tpp_new%plug%edge=plug%edge
 
!  call check_ring('counter ')

 temppx_b=>tpp_new

!  pppl=crpts(1)%temppx%plug
 temppx=>temppx%next
!  pppl=crpts(1)%temppx%plug

end do

! call check_ring('counter ')

! pppl=crpts(1)%temppx%plug

! call check_list_c(7,'down  ')

! Delete the outgoing curve and its plug into the new node cell.
 
! plug=tempp%plug; address=tempp%address
call deletee(cvv(cv_outgoing))

! pppl=crpts(1)%temppx%plug

call deletee(tempp,'after ')

! call check_ring('clock   ')

! address=crpts(1)%temppx%address

! Delete the olde node cell.
call deletee(ndd(nd_nmb))

! call check_ring('counter ')

! call check_list_c(7,'down  ') 

! Recompute the ordinary cell-average in the new node cell.
call evaluate_node_state(nd_nmb_n,or_state)
sum=sum+ndd(nd_nmb_n)%n_cell%or_state-or_state
ndd(nd_nmb_n)%n_cell%or_state=or_state

! Distribute the conservation difference to the critical cells
! located in the cell of the old node cell.
! First decompose the difference in conservation for each curve.
sum=sum/dble(number_of_outgoing)
allocate(dist_cov(number_of_outgoing,wave_number))
do i=1, number_of_outgoing; do j=1, wave_number
 dist_cov(i,j)=0.d0
end do; end do 
allocate(dist_conv(wave_number))
do i=1, wave_number
 dist_conv(i)=0.d0
end do

! pppl=crpts(1)%temppx%plug

do i=1, number_of_outgoing
 crpts_i%temp=>crpts(i)%temp; crpts_i%temppx=>crpts(i)%temppx
 call decompose(crpts_i)
 if(crpts_i%temppx%plug%end_type.eq.'begin ') then
  allocate(xxx(wave_number))
  do j=1, wave_number
   xxx(j)=dist_conv(wave_number+1-j)
  end do
  do j=1, wave_number
   dist_conv(j)=xxx(j)
  end do
  deallocate(xxx)
 end if
 dist_cov(i,:)=dist_conv
end do
! Then distribute it to different waves(curves).

! pppl=crpts(1)%temppx%plug

call distribution

! call check_list_c(7,'down  ')

call change_node_number(nd_nmb_n,nd_nmb)
call renumber_plugs(nd_nmb)
nd_nmb=nd_mem

! call check_ring(1,'counter ')

deallocate(dist_cov); deallocate(dist_conv); deallocate(crpts)

end subroutine awake_curves


subroutine yawn_curve
! This subroutine deals with the case of 'yawn'-status outgoing
! discontinuity curves.

implicit none

type(curve_plug), pointer :: other_tempp, neighbor_tempp, &
                             other_neighbor_tempp, tp_point
type(cv_plug_info) :: plug, ot_plug
integer :: end_node, other_end_node
!real(8), dimension(2) :: posi
type(state) :: l_state, r_state, value
type(phy_info), dimension(wave_number) :: p_new
type(phy_info), dimension(:,:), allocatable :: pp_new
type(geo_info), dimension(:), allocatable :: gg_new
character*3 :: broken, other_broken
integer :: i, j, i0, j0, i_e, j_e, i_oe, j_oe
integer :: cell_number, curve, other_curve, single, exchange, nd_mem
character*8 :: dir
character*12, dimension(wave_number) :: wave_type
logical :: if_first
character*6 :: end_type, other_end_type
character*8 :: insert_order

type(state), dimension(:), allocatable :: xxx

! type(cv_plug_info) :: plugg
! type(adss_plug) :: adsss

other_tempp=>tempp_outgoing%next
plug=tempp_outgoing%plug; ot_plug=other_tempp%plug
curve=plug%cv_nmb; other_curve=ot_plug%cv_nmb

! First, find the two 'yawn'-status curve plugs and the other two
! node cells.
if(cvv(curve)%begin_end.eq.nd_nmb) then
 end_node=cvv(curve)%end_end; end_type='end   '
else
 if(cvv(curve)%end_end.eq.nd_nmb) then
  end_node=cvv(curve)%begin_end; end_type='begin '
 else
  call error_message
 end if
end if
if(cvv(other_curve)%begin_end.eq.nd_nmb) then
 other_end_node=cvv(other_curve)%end_end; other_end_type='end   '
else
 if(cvv(other_curve)%end_end.eq.nd_nmb) then
  other_end_node=cvv(other_curve)%begin_end; other_end_type='begin '
 else
  call error_message
 end if
end if

! The following finds out how many critical cells there will be
! on the connecting curve.
call visit(end_node,curve,end_type,neighbor_tempp)
call visit(other_end_node,other_curve,other_end_type, &
           other_neighbor_tempp)
broken='no '; other_broken='no '; cell_number=1
i0=ndd(nd_nmb)%n_cell%x_idx; j0=ndd(nd_nmb)%n_cell%y_idx
i_e=ndd(end_node)%n_cell%x_idx; j_e=ndd(end_node)%n_cell%y_idx
i_oe=ndd(other_end_node)%n_cell%x_idx
j_oe=ndd(other_end_node)%n_cell%y_idx
if(iabs(i_e-i0)+iabs(j_e-j0).gt.1) then
 broken='yes'; cell_number=cell_number+1
end if
if(iabs(i_oe-i0)+iabs(j_oe-j0).gt.1) then
 other_broken='yes'; cell_number=cell_number+1
end if
allocate(pp_new(wave_number,cell_number))
allocate(gg_new(cell_number))

! The following finds out geometrical information for the new 
! critical cells.
i=1

if(broken.eq.'yes') then
! Set the critical cell that fixes the broken between the end
! node cell and the presetn node cell.
 gg_new(i)%edge(1)=neighbor_edge(neighbor_tempp%plug%edge)
 gg_new(i)%edge(2)=neighbor_edge(plug%edge)
 call find_points_from_edges(gg_new(i))
 call find_type_from_edges(gg_new(i))
 select case(gg_new(i)%edge(2))
  case(1)
   gg_new(i)%x_idx=i0; gg_new(i)%y_idx=j0+1
  case(2)
   gg_new(i)%x_idx=i0-1; gg_new(i)%y_idx=j0
  case(3)
   gg_new(i)%x_idx=i0; gg_new(i)%y_idx=j0-1
  case(4)
   gg_new(i)%x_idx=i0+1; gg_new(i)%y_idx=j0
 end select
 call find_single(gg_new(i),single)
 select case(single)
  case(1); gg_new(i)%dis_posi=-0.4999999
  case(2)
   if(gg_new(i)%edge(1).eq.1) then
    gg_new(i)%dis_posi(1)=0.4999999
    gg_new(i)%dis_posi(2)=-0.4999999
   else
    gg_new(i)%dis_posi(1)=-0.4999999
    gg_new(i)%dis_posi(2)=0.4999999
   end if
  case(3); gg_new(i)%dis_posi=0.4999999
  case(4)
   if(gg_new(i)%edge(1).eq.4) then
    gg_new(i)%dis_posi(1)=0.4999999
    gg_new(i)%dis_posi(2)=-0.4999999
   else
    gg_new(i)%dis_posi(1)=-0.4999999
    gg_new(i)%dis_posi(2)=0.4999999
   end if
 end select
 i=i+1
end if

! Set the critical cell locating in the grid cell of the present
! node cell.
gg_new(i)%x_idx=i0; gg_new(i)%y_idx=j0
gg_new(i)%edge(1)=plug%edge; gg_new(i)%edge(2)=ot_plug%edge
if(broken.eq.'yes') then
 gg_new(i)%dis_posi(1)=gg_new(i-1)%dis_posi(2)
else 
 gg_new(i)%dis_posi(1)=0.0d0
end if
gg_new(i)%dis_posi(2)=0.d0
call find_type_from_edges(gg_new(i))
call find_points_from_edges(gg_new(i))
i=i+1

if(other_broken.eq.'yes') then
! Set the critical cell that fixes the broken between the other
! end node and the present node.
 gg_new(i)%edge(1)=neighbor_edge(ot_plug%edge)
 gg_new(i)%edge(2)=neighbor_edge(other_neighbor_tempp%plug%edge)
 call find_points_from_edges(gg_new(i))
 call find_type_from_edges(gg_new(i))
 select case(gg_new(i)%edge(1))
  case(1)
   gg_new(i)%x_idx=i0; gg_new(i)%y_idx=j0+1
  case(2)
   gg_new(i)%x_idx=i0-1; gg_new(i)%y_idx=j0
  case(3)
   gg_new(i)%x_idx=i0; gg_new(i)%y_idx=j0-1
  case(4)
   gg_new(i)%x_idx=i0+1; gg_new(i)%y_idx=j0
 end select
 call find_single(gg_new(i),single)
 select case(single)
  case(1); gg_new(i)%dis_posi=-0.4999999
  case(2)
   if(gg_new(i)%edge(1).eq.1) then
    gg_new(i)%dis_posi(1)=0.4999999
    gg_new(i)%dis_posi(2)=-0.4999999
   else
    gg_new(i)%dis_posi(1)=-0.4999999
    gg_new(i)%dis_posi(2)=0.4999999
   end if
  case(3); gg_new(i)%dis_posi=0.4999999
  case(4)
   if(gg_new(i)%edge(1).eq.4) then
    gg_new(i)%dis_posi(1)=0.4999999
    gg_new(i)%dis_posi(2)=-0.4999999
   else
    gg_new(i)%dis_posi(1)=-0.4999999
    gg_new(i)%dis_posi(2)=0.4999999
   end if
 end select
 gg_new(i-1)%dis_posi(2)=gg_new(i)%dis_posi(1)
end if

! The following sets the physical information for the new 
! critical cell.

! Set left and right states and the ring direction of the
! corresponding plugs.
i=1
if(broken.eq.'yes') then
 l_state=plug%u_behind(gg_new(i)%x_idx-x_nd,gg_new(i)%y_idx-y_nd)
 r_state=plug%u_front(gg_new(i)%x_idx-x_nd,gg_new(i)%y_idx-y_nd)
 dir='clock   '; if_first=.true.
 call physical_separation(l_state,r_state,gg_new(i),dir,p_new, &
      wave_type,if_first,number_of_outgoing)
 do	j=1, wave_number; pp_new(j,i)=p_new(i); end do
 i=i+1
end if

l_state=plug%u_behind(0,0)
r_state=plug%u_front(0,0)
dir='clock   '
if(broken.eq.'yes') then
 if_first=.false.
else 
 if_first=.true.
endif
call physical_separation(l_state,r_state,gg_new(i),dir,p_new, &
     wave_type,if_first,number_of_outgoing)
do j=1, wave_number; pp_new(j,i)=p_new(j); end do
i=i+1

if(other_broken.eq.'yes') then
 l_state=plug%u_behind(gg_new(i)%x_idx-x_nd,gg_new(i)%y_idx-y_nd)
 r_state=plug%u_front(gg_new(i)%x_idx-x_nd,gg_new(i)%y_idx-y_nd)
 dir='clock   '; if_first=.false.
 call physical_separation(l_state,r_state,gg_new(i),dir,p_new, &
      wave_type,if_first,number_of_outgoing)
 do	j=1, wave_number; pp_new(j,i)=p_new(j); end do
 i=i+1
end if

! call check_list_c(2,'down  ') 
! call check_ring(1,'counter ')
 
! The following updates the numerical solution.
! First delete the present node cell and let the two new node cells
! be connected by 'curve'.
call deletee(ndd(nd_nmb))
nd_mem=nd_nmb
exchange=nd_nmb; nd_nmb=end_node; end_node=exchange
nd_nmb_n=other_end_node 
neighbor_tempp%plug%end_type='begin '
other_neighbor_tempp%plug%end_type='end   '
other_neighbor_tempp%plug%cv_nmb=curve 
cvv(curve)%begin_end=nd_nmb; cvv(curve)%end_end=nd_nmb_n
call deletee(cvv(other_curve))

tempp=>neighbor_tempp

call yawn_update(gg_new,pp_new,cell_number,insert_order)

! call check_ring('counter ')
! call check_list_c(5,'down  ')

allocate(crpts(wave_number))
call evaluation_in_cell(i0,j0,value)
sum=sum-value
tp_point=>tempp
do i=1,number_of_outgoing
 select case(insert_order)
  case('counter '); j=number_of_outgoing+1-i
  case('clock   '); j=i
  case default; call error_message
 end select
 crpts(j)%temppx=>tp_point
 if(broken.eq.'yes') then
  crpts(j)%temp=>cvv(tp_point%plug%cv_nmb)%begin%next%next 
 else
  crpts(j)%temp=>cvv(tp_point%plug%cv_nmb)%begin%next
 end if
 select case(insert_order)
  case('counter '); tp_point=>tp_point%previous
  case('clock   '); tp_point=>tp_point%next
  case default; call error_message
 end select     
end do

! call check_list_c(5,'down  ')

! plugg=crpts(1)%temppx%plug; adsss=crpts(1)%temppx%address
sum=sum/dble(number_of_outgoing)
allocate(dist_cov(number_of_outgoing,wave_number))
do i=1, number_of_outgoing; do j=1, wave_number
 dist_cov(i,j)=0.d0
end do; end do 
allocate(dist_conv(wave_number))
do i=1, wave_number
 dist_conv(i)=0.d0
end do

! pppl=crpts(1)%temppx%plug

do i=1, number_of_outgoing
 crpts_i%temp=>crpts(i)%temp; crpts_i%temppx=>crpts(i)%temppx
 call decompose(crpts_i)
 if(crpts_i%temppx%plug%end_type.eq.'end   ') then
  allocate(xxx(wave_number))
  do j=1, wave_number
   xxx(j)=dist_conv(wave_number+1-j)
  end do
  do j=1, wave_number
   dist_conv(j)=xxx(j)
  end do
  deallocate(xxx)
 end if
 dist_cov(i,:)=dist_conv
end do
! Then distribute it to different waves(curves).

! pppl=crpts(1)%temppx%plug

! call check_list_c(5,'down  ')

call distribution

! call check_list_c(5,'down  ')
call change_node_number(nd_nmb,end_node)
call renumber_plugs(nd_nmb)
nd_nmb=nd_mem

deallocate(gg_new); deallocate(pp_new)
deallocate(dist_cov); deallocate(dist_conv); deallocate(crpts)

end subroutine yawn_curve


subroutine insert_n_update(if_in_old_cell)
! This subroutine carries out the job of insertion new critical cell
! and update the numerical solution in the vicinity.

implicit none
character*3, intent(in) :: if_in_old_cell

! type(geo_info) :: g_cell

! g_cell=temp%g_cell

allocate(t_new); t_new%g_cell=g_new; t_new%p_cell=p_new
select case(plug%end_type)
 case('begin '); call insert(temp,t_new,'before')
 case('end   '); call insert(temp,t_new,'after ')
end select

temp=>t_new

! call check_list_c(3,'up    ')

! Finally, update the solution in smooth, update the grid-map and
! preserve conservation.
call update_near_node(t_new,tempp)

! The conservation-preserving operation is carried out only in
! the critical cell between the head of discontinuity curve and 
! and the old node cell and the between the old node cell and 
! the new node cell
select case(if_in_old_cell) 
 case('yes')
  sum=sum-p_new%or_state
  if(t_new%l_stk.ne.adss_info(0,0)) sum=sum+p_new%l_state
  if(t_new%r_stk.ne.adss_info(0,0)) sum=sum+p_new%r_state
  if(t_new%l_stk.ne.adss_info(0,0).and.t_new%r_stk.ne.adss_info(0,0)) call error_message
 case('no ')
  call find_single(g_new,single)
  select case(g_new%point(single))
   case('left  '); dif_cov=p_new%or_state-p_new%r_state
   case('right '); dif_cov=p_new%or_state-p_new%l_state
   case default; print*, 'Something is wrong!!!'; pause
  end select
  sum=sum-dif_cov

 case default; print*, 'Something is wrong'; pause
end select

end subroutine insert_n_update


subroutine set_start
! This subroutine sets the start information

implicit none

ed_posi(1)=ndd(nd_nmb_n)%n_cell%x_posi+x_nd_new-g_new%x_idx
ed_posi(2)=ndd(nd_nmb_n)%n_cell%y_posi+y_nd_new-g_new%y_idx
select case(plug%end_type)
 case('begin ')
  call pick_point(g_cell,1,st_posi); st_edge=g_cell%edge(1)
 case('end   ')
  call pick_point(g_cell,2,st_posi); st_edge=g_cell%edge(2)
 case default; print*, 'Something is wrong!!!'; pause
end select
st_edge=neighbor_edge(st_edge)
st_posi(1)=st_posi(1)+dfloat(g_cell%x_idx-g_new%x_idx)
st_posi(2)=st_posi(2)+dfloat(g_cell%y_idx-g_new%y_idx)
ed_edge=error_index

end subroutine set_start


subroutine decompose(crpts_i)
! This subroutine decomposes the difference in conservation for
! a give curve (in a critical cell).

implicit none
!integer, intent(in) :: i
type(cr_pointer), intent(in) :: crpts_i

type(geo_info) :: gd_cell
type(phy_info) :: pd_cell
real(8), dimension(2) :: normal, pt1, pt2
type(state) :: l_state, r_state
type(cv_plug_info) :: pdlug
integer :: wave_nb

gd_cell=crpts_i%temp%g_cell; pd_cell=crpts_i%temp%p_cell

! Compute the normal of the discontinuity curve in the critical
! cell first.

call pick_point(gd_cell,1,pt1); call pick_point(gd_cell,2,pt2)
normal=normal_of_line(pt1,pt2)

! Determine the direction of normal and then take the physical 
! states on the two sides accordingly.

pdlug=crpts_i%temppx%plug; wave_nb=pd_cell%wv_nb

!select case(pdlug%end_type)
! case('begin ')
  l_state=pd_cell%l_state; r_state=pd_cell%r_state
! case('end   ')
!  l_state=pd_cell%r_state; r_state=pd_cell%l_state
!  normal=-normal; wave_nb=wave_number+1-wave_nb
! case default; print*, 'Something is wrong!!!'; pause
!end select

! Decompose the conservation difference
call decompose_conservation(l_state,r_state,normal,sum,wave_nb,dist_conv)

end subroutine decompose


subroutine distribution
! This subroutine distributes the conservation difference among
! the outgoing discontinuity curves.

implicit none
type(state), dimension(wave_number) :: conv
integer :: i, j, ja1, ja2, jb1, jb2, jj
type(cv_plug_info) :: plug_a, plug_b, plug
type(phy_info), dimension(:), allocatable :: pp_cell
type(state) :: a, b
real(8) :: coef_a, coef_b

! type(geo_info) :: g_cell

allocate(pp_cell(number_of_outgoing))
do i=1, number_of_outgoing
 pp_cell(i)=crpts(i)%temp%p_cell
end do

do i=1,wave_number
 conv(i)=0.d0
 do j=1,number_of_outgoing
  conv(i)=conv(i)+dist_cov(j,i)
 end do
end do

if(wave_number.eq.1) then
 if(number_of_outgoing.gt.1) call error_message
 pp_cell(1)%or_state=pp_cell(1)%or_state+conv(1)
 crpts(1)%temp%p_cell=pp_cell(1)
 return
end if
 
do i=1,number_of_outgoing
 a=0.0d0; b=0.0d0
 plug=crpts(i)%temppx%plug
 if(plug%end_type.eq.'begin ') then
  jj=wave_number-pp_cell(i)%wv_nb+1
 else
  jj=pp_cell(i)%wv_nb
 end if

 if(i.eq.1) then
  ja1=1; call clean_up_plug(plug_a); coef_a=1.0d0
 else
  plug_a=crpts(i-1)%temppx%plug 
  if(plug_a%end_type.eq.'begin ') then
   ja1=wave_number-pp_cell(i-1)%wv_nb+2
  else
   ja1=pp_cell(i-1)%wv_nb+1
  end if
!  ja1=plug_a%outgoing_wave_number+1; coef_a=0.5d0
  coef_a=0.5d0
 end if
 if(plug%end_type.eq.'begin ') then
  ja2=wave_number-pp_cell(i)%wv_nb
 else
  ja2=pp_cell(i)%wv_nb-1
 end if
 !ja2=plug%outgoing_wave_number-1
 do j=ja1,ja2
  a=a+conv(j)
 end do

 if(i.eq.number_of_outgoing) then
  jb2=wave_number; call clean_up_plug(plug_b); coef_b=1.0d0
 else
  plug_b=crpts(i+1)%temppx%plug
  if(plug_b%end_type.eq.'begin') then
   jb2=wave_number-pp_cell(i+1)%wv_nb
  else
   jb2=pp_cell(i+1)%wv_nb-1
  end if
!  jb2=plug_b%outgoing_wave_number-1; coef_b=0.5d0
 end if
 if(plug%end_type.eq.'begin ') then
  jb1=wave_number-pp_cell(i)%wv_nb+2
 else
  jb1=pp_cell(i)%wv_nb+1
 end if
! jb1=plug%outgoing_wave_number+1
 do j=jb1,jb2
  b=b+conv(j)
 end do

 select case(crpts(i)%temppx%plug%end_type) 
  case('begin ')
   pp_cell(i)%l_state=pp_cell(i)%l_state+coef_b*b
   pp_cell(i)%r_state=pp_cell(i)%r_state+coef_a*a
  case('end   ')
   pp_cell(i)%l_state=pp_cell(i)%l_state+coef_a*a
   pp_cell(i)%r_state=pp_cell(i)%r_state+coef_b*b
  case default; print*, 'Something is wrong!!!'; pause
 end select
 pp_cell(i)%or_state=pp_cell(i)%or_state+coef_a*a+coef_b*b+conv(jj)
end do

do i=1, number_of_outgoing

! g_cell=crpts(i)%temp%g_cell

 crpts(i)%temp%p_cell=pp_cell(i)
end do

deallocate(pp_cell)

end subroutine distribution


end module fix_triple_nodes
