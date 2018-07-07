module tools_4_connecting_n_triple
! This module provides tools necessary for dealing with fixing 
! broken connections between discontinuity curves and the node
! cell they are supposed to plug-in and triple node cells.

use variables_4_node_treat
! 'bs_var.f90'

implicit none

public  new_geo_cell_n_area, physical_states, & 
        new_geo_cell_4_yawn, yawn_update
private find_cuts


contains


subroutine new_geo_cell_n_area(st_posi,st_edge,ed_posi,ed_edge, &
           end_type,g_new,area,edd_edge)
! Given start discontinuity position (edge plus position), end
! edge and end discontinuity position, this subroutine finds 
! the geometrical information for a new critical cell that starts 
! discontinuity position, ends at the end edge, and is cut throuth
! by the line segment connecting the start and end discontinuity
! positions and also gives the area cut by the segment.

implicit none
real(8), dimension(2), intent(in) :: st_posi,ed_posi
integer, intent(in) :: st_edge
integer, intent(inout) :: ed_edge
integer, intent(out) :: edd_edge
character*6, intent(in) :: end_type
type(geo_info), intent(inout) ::g_new
real(8), intent(out) :: area

type(polygon) :: p1, p2, p3
real(8), dimension(2) :: pt1, pt2
real(8), dimension(2) :: cut
logical :: manipulated

manipulated=.false.
! The second polygon is taken to be the new critical cell.
p2%x(1)=-0.5d0; p2%y(1)=-0.5d0; p2%x(2)=-0.5d0; p2%y(2)=0.5d0
p2%x(3)=0.5d0; p2%y(3)=0.5d0; p2%x(4)=0.5d0; p2%y(4)=-0.5d0
p2%size=4

! The first polygon is taken to be the discontinuity segment.
p1%size=2
select case(end_type)
 case('begin ')
  p1%x(2)=st_posi(1)
  p1%y(2)=st_posi(2)
  p1%x(1)=ed_posi(1)
  p1%y(1)=ed_posi(2)   
 case('end   ')
  p1%x(1)=st_posi(1)
  p1%y(1)=st_posi(2)
  p1%x(2)=ed_posi(1)
  p1%y(2)=ed_posi(2)   
 case default; print*, 'Something wrong!!!'; pause
end select
pt1(1)=p1%x(1); pt1(2)=p1%y(1); pt2(1)=p1%x(2); pt2(2)=p1%y(2)

! The following finds out the requested information.
call polyin(p1,p2,p3,'open  ',area)
call find_cuts(pt1,pt2,p2,st_edge,cut,ed_edge,manipulated)
if(manipulated) then
 select case(st_edge)
  case(1); edd_edge=3
  case(2); edd_edge=4
  case(3); edd_edge=1
  case(4); edd_edge=2
 end select
end if

! Compute geometrical information of the new critical cell.
select case(end_type)
 case('begin ')
  g_new%edge(1)=ed_edge; g_new%edge(2)=st_edge
  select case(g_new%edge(1))
   case(1,3); g_new%dis_posi(1)=cut(1)
   case(2,4); g_new%dis_posi(1)=cut(2)
  end select
!  if(dabs(g_new%dis_posi(1)).gt.0.5d0) then
!   g_new%dis_posi(1)=dsign(0.4999999d0,g_new%dis_posi(1))
!  end if
  select case(g_new%edge(2))
   case(1,3); g_new%dis_posi(2)=st_posi(1)
   case(2,4); g_new%dis_posi(2)=st_posi(2)
  end select
 case('end   ')
  g_new%edge(1)=st_edge; g_new%edge(2)=ed_edge
  select case(g_new%edge(1))
   case(1,3); g_new%dis_posi(1)=st_posi(1)
   case(2,4); g_new%dis_posi(1)=st_posi(2)
  end select
  select case(g_new%edge(2))
   case(1,3); g_new%dis_posi(2)=cut(1)
   case(2,4); g_new%dis_posi(2)=cut(2)
  end select
!  if(dabs(g_new%dis_posi(2)).gt.0.5d0) then
!   g_new%dis_posi(2)=dsign(0.4999999d0,g_new%dis_posi(2))
!  end if
end select

call find_points_from_edges(g_new)
call find_type_from_edges(g_new)

end subroutine new_geo_cell_n_area


subroutine find_cuts(pt1,pt2,poly,st_edge,cut,ed_edge,manipulated)
! This subroutine finds the two cut-points of a line segment with
! a polygon.

implicit none
real(8), dimension(2), intent(in) :: pt1, pt2
! The two end-points of the line segment.
type(polygon), intent(in) :: poly
! The polygon under concern.
real(8), dimension(2), intent(out) :: cut
! The two cuts cut by the line segment from the polygon.
integer, intent(in) :: st_edge
integer, intent(inout) :: ed_edge
! The edges the cut is on.
logical, intent(out) :: manipulated
character*10 :: error_msg

integer :: i, ll
real(8), dimension(2,2) :: pts1, pts2
real(8), dimension(2) :: pt, ptt1, ptt2
real(8), dimension(2) :: x_dif, y_dif
integer :: start_edge, end_edge, side1, side2
real(8) :: a, b, c, dis1, dis2

real(8), dimension(2,2) :: ccut
integer, dimension(2) :: eend_edge
character*3, dimension(2) :: if_inside
real(8) :: distance_1, distance_2

manipulated=.false.
pts1(1,:)=pt1; pts1(2,:)=pt2

start_edge=5-st_edge; end_edge=5-ed_edge
! The vortices of the ploygon are numbered in clockwise way while
! its edges are numbered in counter-clockwise way, thus such a
! transform of numbers of edges is needed.

select case(ed_edge)

 case (1,2,3,4)
! If the end-edge is defined.
  pts2(1,1)=poly%x(end_edge); pts2(1,2)=poly%y(end_edge)

  if(end_edge.lt.poly%size) then
   pts2(2,1)=poly%x(end_edge+1); pts2(2,2)=poly%y(end_edge+1)
  else
   pts2(2,1)=poly%x(1); pts2(2,2)=poly%y(1)
  end if
  ptt1=pts1(1,:); ptt2=pts1(2,:)
  call line_pass_points(ptt1,ptt2,a,b,c)
  pt=pts2(1,:); dis1=distance_to_line(a,b,c,pt)
  pt=pts2(2,:); dis2=distance_to_line(a,b,c,pt)
  if(dis1*dis2.ge.0.d0) then; side2=1; else; side2=-1; end if
  ptt1=pts2(1,:); ptt2=pts2(2,:)
  call line_pass_points(ptt1,ptt2,a,b,c)
  pt=pts1(1,:); dis1=distance_to_line(a,b,c,pt)
  pt=pts1(2,:); dis2=distance_to_line(a,b,c,pt)
  if(dis1*dis2.ge.0.d0) then; side1=1; else; side1=-1; end if
  if(side1.lt.0.or.side2.lt.0) then
   call ints_of_lines(pts1,pts2,pt,error_msg)
   do i=1,2
    if(pt(i).lt.-0.5d0) pt(i)=-0.4999999
    if(pt(i).gt.0.5d0) pt(i)=0.4999999
   end do
   cut(1)=pt(1); cut(2)=pt(2)
  else
   pt=pts1(1,:)
   dis1=length_of_segment(ptt1,pt)
   dis2=length_of_segment(ptt2,pt)
   manipulated=.true.
   if(dis1.lt.dis2) then
    cut=ptt1
!    if(end_edge.lt.poly%size) then
!     edd_edge=end_edge+1
!    else
!     edd_edge=1
!    end if	   	     
   else
    cut=ptt2
!    if(end_edge.gt.1) then
!     edd_edge=end_edge-1
!    else
!     edd_edge=poly%size
!    end if	   	     
   end if
   select case(ed_edge)
    case(1,3); cut(1)=0.4999999*dabs(cut(1))/cut(1)
    case(2,4); cut(2)=0.4999999*dabs(cut(2))/cut(2)
   end select
  end if 

case default
! If the end-edge is not defined.
  ccut=error_data; eend_edge=error_index; ll=1
  do i=1, poly%size
   if(i.eq.start_edge) cycle
! The start-edge should be excluded.
   pts2(1,1)=poly%x(i); pts2(1,2)=poly%y(i)
   if(i.lt.poly%size) then
    pts2(2,1)=poly%x(i+1); pts2(2,2)=poly%y(i+1)
   else
    pts2(2,1)=poly%x(1); pts2(2,2)=poly%y(1)
   end if
   call ints_of_lines(pts1,pts2,pt,error_msg)
   x_dif(1)=pt(1)-pts1(1,1); y_dif(1)=pt(2)-pts1(1,2)
   x_dif(2)=pt(1)-pts1(2,1); y_dif(2)=pt(2)-pts1(2,2)
   if(x_dif(1)*x_dif(2).lt.0.d0.or.y_dif(1)*y_dif(2).lt.0.d0) then
    if(ll.gt.2) then
     print*, 'Something wrong!!!'; pause
    end if
    ccut(ll,1)=pt(1); ccut(ll,2)=pt(2); eend_edge(ll)=i
    ll=ll+1
   end if
  end do
  ll=ll-1

  select case(ll)
   case(1)
    pts2(1,1)=poly%x(eend_edge(1)); pts2(1,2)=poly%y(eend_edge(1))
    if(eend_edge(1).lt.poly%size) then
     pts2(2,1)=poly%x(eend_edge(1)+1)
     pts2(2,2)=poly%y(eend_edge(1)+1)
    else
     pts2(2,1)=poly%x(1); pts2(2,2)=poly%y(1)
    end if
    x_dif(1)=ccut(1,1)-pts2(1,1); y_dif(1)=ccut(1,2)-pts2(1,2)
    x_dif(2)=ccut(1,1)-pts2(2,1); y_dif(2)=ccut(1,2)-pts2(2,2)
    if(x_dif(1)*x_dif(2).gt.0.d0.or.y_dif(1)*y_dif(2).gt.0.d0) then
     distance_1=dsqrt(x_dif(1)*x_dif(1)+y_dif(1)*y_dif(1))
     distance_2=dsqrt(x_dif(2)*x_dif(2)+y_dif(2)*y_dif(2))
     if(distance_1.lt.distance_2) then
	  ccut(1,:)=pts2(1,:)+0.0000002*(pts2(2,:)-pts2(1,:))
     else
	  ccut(1,:)=pts2(2,:)+0.0000002*(pts2(1,:)-pts2(2,:))
     end if
	end if
	cut=ccut(1,:); end_edge=eend_edge(1)

   case(2)
    do i=1,2
     pts2(1,1)=poly%x(eend_edge(i)); pts2(1,2)=poly%y(eend_edge(i))
     if(eend_edge(i).lt.poly%size) then
      pts2(2,1)=poly%x(eend_edge(i)+1)
      pts2(2,2)=poly%y(eend_edge(i)+1)
     else
      pts2(2,1)=poly%x(1); pts2(2,2)=poly%y(1)
     end if
     x_dif(1)=ccut(i,1)-pts2(1,1); y_dif(1)=ccut(i,2)-pts2(1,2)
     x_dif(2)=ccut(i,1)-pts2(2,1); y_dif(2)=ccut(i,2)-pts2(2,2)
     if(x_dif(1)*x_dif(2).lt.0.d0.or.y_dif(1)*y_dif(2).lt.0.d0) then
      if_inside(i)='yes'
     else
      if_inside(i)='no '
     end if
    end do
    if(if_inside(1).eq.'yes'.and.if_inside(2).eq.'yes') call error_message
    if(if_inside(1).eq.'no '.and.if_inside(2).eq.'no ') call error_message
    if(if_inside(1).eq.'yes') then
     cut(1)=ccut(1,1); cut(2)=ccut(1,2); end_edge=eend_edge(1)
    else
     cut(1)=ccut(2,1); cut(2)=ccut(2,2); end_edge=eend_edge(2)
    end if

  end select

end select 

ed_edge=5-end_edge

end subroutine find_cuts


subroutine physical_states(plug,g_cell,area,p_cell)
! This subroutine computes the physical states for the new 
! critical cell.

implicit none
type(cv_plug_info), intent(in) :: plug
type(geo_info), intent(in) :: g_cell
real(8), intent(in) :: area
type(phy_info), intent(out) :: p_cell

select case(plug%end_type)
 case('begin ')
  p_cell%l_state=plug%u_front(g_cell%x_idx-x_nd,g_cell%y_idx-y_nd)
  p_cell%r_state=plug%u_behind(g_cell%x_idx-x_nd,g_cell%y_idx-y_nd)
 case('end   ')
  p_cell%l_state=plug%u_behind(g_cell%x_idx-x_nd,g_cell%y_idx-y_nd)
  p_cell%r_state=plug%u_front(g_cell%x_idx-x_nd,g_cell%y_idx-y_nd)
end select
p_cell%or_state=area*p_cell%r_state+(1.d0-area)*p_cell%l_state

end subroutine physical_states


subroutine new_geo_cell_4_yawn(st_edge,ed_edge,end_type,posi,g_new)
! Given start edge-plus and end-edge, this subroutine finds 
! the geometrical information for a new critical cell that is one 
! of the first two critical cell s for a 'yawn' status critical
! cell.

implicit none
integer, intent(in) :: st_edge
integer, intent(inout) :: ed_edge
character*6, intent(in) :: end_type
real(8), dimension(2) :: posi
type(geo_info), intent(inout) ::g_new

integer :: num_b, num_e

! Set numbers of edges
select case(end_type)
 case('begin '); num_b=2; num_e=1
 case('end   '); num_b=1; num_e=2
 case default; print*, 'Something is wrong!!!'; pause
end select
! Set edges for the new critical cell.
g_new%edge(num_b)=ed_edge
g_new%edge(num_e)=st_edge
! Set discontinuity positions for the new critical cell.
select case(g_new%edge(num_b))
 case(1)
  g_new%dis_posi(num_e)=-posi(2)
  select case(x_nd-g_new%x_idx)
   case(1)
    g_new%dis_posi(num_b)=posi(1)
   case(-1)
    g_new%dis_posi(num_b)=-posi(1)
   case default; print*, 'Something is wrong!!!'; pause
  end select
 case(2)
  g_new%dis_posi(num_e)=posi(2)
  select case(y_nd-g_new%y_idx)
   case(1)
    g_new%dis_posi(num_b)=posi(1)
   case(-1)
    g_new%dis_posi(num_b)=-posi(1)
   case default; print*, 'Something is wrong!!!'; pause
  end select
 case(3)
  g_new%dis_posi(num_e)=posi(2)
  select case(x_nd-g_new%x_idx)
   case(1)
    g_new%dis_posi(num_b)=posi(1)
   case(-1)
    g_new%dis_posi(num_b)=-posi(1)
   case default; print*, 'Something is wrong!!!'; pause
  end select
 case(4)
  g_new%dis_posi(num_e)=-posi(2)
  select case(y_nd-g_new%y_idx)
   case(1)
    g_new%dis_posi(num_b)=posi(1)
   case(-1)
    g_new%dis_posi(num_b)=-posi(1)
   case default; print*, 'Something is wrong!!!'; pause
  end select      
end select
! Set the type and points for the new critical cell.
call find_points_from_edges(g_new)
call find_type_from_edges(g_new)

end subroutine new_geo_cell_4_yawn


subroutine physical_separation(l_state,r_state,g_cell,dir, &
           p_cell,wave_type,if_first,dis_number) 
! Given geometrical information of a critical cell and a left and
! a right states this subroutine solves the Riemann problem to
! separate the left and right states and forms corresponding 
! physical critical cells. This subroutine is used in producing
! critical cells for a 'yawn' status discontinuity
! curves

implicit none
type(state), intent(in) :: l_state, r_state
type(geo_info), intent(in) :: g_cell
character*8, intent(in) :: dir
! Indicate if the plugs of the three discontinuity curves are in
! clockwise or counter-clockwise order.
type(phy_info), dimension(wave_number), intent(out) :: p_cell
logical, intent(in) :: if_first
character*12, dimension(wave_number), intent(inout) :: wave_type
integer, intent(out) :: dis_number

type(state) :: c1, c2
real(8), dimension(2) :: pt1, pt2, normal
integer :: i, ii, iii, wave_before, wave_after
character*12 :: wv_type
real(8) :: sp, area_l, area_r, x_posi, y_posi

if(if_first) wave_type='rarefaction '
wave_before=0; wave_after=0; dis_number=0; iii=0
do ii=1, wave_number
 call clean_up_p_cell(p_cell(ii))
 select case(dir)
  case('clock   '); i=wave_number+1-ii
  case('counter '); i=ii 
 end select
 call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2)
 normal=normal_of_line(pt1,pt2)
 x_posi=dfloat(g_cell%x_idx)*h; y_posi=dfloat(g_cell%y_idx)*h
 call riemann(l_state,r_state,normal,i,x_posi,y_posi,c1,c2,sp,wv_type)
 if(if_first) wave_type(ii)=wv_type
 if(wave_type(ii).eq.'shock       '.or.wave_type(ii).eq.'contact     ') then
  iii=iii+1
  if(wave_before.gt.0) then
   p_cell(iii)%l_state=0.5d0*(c1+p_cell(wave_before)%r_state)
   p_cell(wave_before)%r_state=p_cell(iii)%l_state
  else
   p_cell(iii)%l_state=c1
  end if
  if(wave_after.gt.0) then
   p_cell(iii)%r_state=0.5d0*(c2+p_cell(wave_after)%l_state)
   p_cell(wave_after)%l_state=p_cell(iii)%r_state
  else
   p_cell(iii)%r_state=c2
  end if
  p_cell(iii)%wv_nb=i
  wave_before=0; wave_after=0
  select case(dir)
   case('clock   '); wave_after=iii
   case('counter '); wave_before=iii
  end select
    
  area_l=side_area(g_cell,'left  ')
  area_r=side_area(g_cell,'right ')
  p_cell(iii)%or_state= &
  area_l*p_cell(iii)%l_state+area_r*p_cell(iii)%r_state
  dis_number=dis_number+1
 end if
end do   

! p_cell%l_state=l_state; p_cell%r_state=r_state

end subroutine physical_separation


subroutine yawn_update(gg_cell,pp_cell,num,insert_order)
! In fixxing broken connections of 'yawn' status discontinuity
! curves and handling triple node cells with 'yawn' status 
! 'outgoing' discontinuity curves, the 'yawn' status curves
! will be separated into several discontinuity curves. Then 
! the 'yawn' curve will be deleted and the new curves will be
! pluged in. This subroutine is to do this delete-plug job.

implicit none
integer,intent(in) :: num
! The number of the critical cells on each new curve.
type(geo_info), dimension(num), intent(in) :: gg_cell
! The geometrical critical cell.
type(phy_info), dimension(wave_number,num), intent(in) :: pp_cell
! The physical critical cells.
character*8, intent(out) :: insert_order

type(node_cell) :: other_node
type(critical_cell), pointer :: temp, t_new
type(curve_plug), pointer :: tp_new, other_tempp, other_tp_new, &
  							 tp_insert, other_tp_insert
type(cv_plug_info) :: plug, pl_new, other_plug, other_new_plug
integer :: i, ii, new_curve, head_mark
logical :: head_switch
character*8 :: curves_order
integer, dimension(:), allocatable :: wave_nb

! type(cv_plug_info) :: plugg
! type(adss_plug) :: adsss

!   plugg=tempp%plug; adsss=tempp%address

! Delete the old plug first.

plug=tempp%plug
! The following finds the other node cell connecting to the curve 
! and the corresponding curve plug.
if(cvv(plug%cv_nmb)%begin_end.eq.nd_nmb) then
 other_node=ndd(cvv(plug%cv_nmb)%end_end)
 call visit(cvv(plug%cv_nmb)%end_end,plug%cv_nmb,'end   ', &
            other_tempp)
end if
if(cvv(plug%cv_nmb)%end_end.eq.nd_nmb) then
 other_node=ndd(cvv(plug%cv_nmb)%begin_end)
 call visit(cvv(plug%cv_nmb)%begin_end,plug%cv_nmb,'begin ', &
            other_tempp)
end if
other_plug=other_tempp%plug

select case(plug%end_type)
 case('begin ') 
  call deletee(tempp,'after ')
  call deletee(other_tempp,'before')
 case('end   ')
  call deletee(tempp,'before')
  call deletee(other_tempp,'after ')
 case default; call error_message
end select
tp_insert=>tempp; other_tp_insert=>other_tempp

! The following decides the the insert order.
curves_order='increase'
if(wave_number.gt.1) then
 allocate(wave_nb(wave_number)); wave_nb=error_index
 do ii=1, wave_number
  if(pp_cell(ii,1)%wv_nb.gt.0) then
   wave_nb(ii)=pp_cell(ii,1)%wv_nb
  end if
 end do
 if(wave_nb(1).gt.0.and.wave_nb(2).gt.0) then
  if(wave_nb(1).gt.wave_nb(2)) curves_order='decrease'
 end if
 deallocate(wave_nb)
end if
 
! The following decides the insert order.
insert_order='counter '
if(plug%end_type.eq.'begin '.and.curves_order.eq.'increase') then
 insert_order='clock   '
end if
if(plug%end_type.eq.'end   '.and.curves_order.eq.'decrease') then
 insert_order='clock   '
end if

! Creat the new discontinuity curves and then plug them in the 
! node cell.
do ii=1,wave_number
 if(pp_cell(ii,1)%wv_nb.gt.0) then
! The discontinuity is produced on the condition that there is
! a corresponding discontinuity resulted in the Riemann separation.

  do i=1, cvn
   if(cvv(i)%status.eq.'asleep') exit
  end do
  new_curve=i
  cvv(new_curve)%status='awake '
  cvv(new_curve)%cv_type=cvv(plug%cv_nmb)%cv_type
  cvv(new_curve)%begin_end=cvv(plug%cv_nmb)%begin_end
  cvv(new_curve)%end_end=cvv(plug%cv_nmb)%end_end
  cvv(new_curve)%wave=pp_cell(ii,1)%wv_nb; cvv(new_curve)%total=0
  call creat_cv(new_curve)
! The above produces the curve and the following inserts the
! critical cell into the curve.
  if(plug%end_type.eq.'begin ') then
   temp=>cvv(new_curve)%begin
  else
   temp=>cvv(new_curve)%eend
  end if
  do i=1, num
   allocate(t_new)
   call clean_up_address(t_new%l_stk); call clean_up_address(t_new%r_stk)
   t_new%g_cell=gg_cell(i); t_new%p_cell=pp_cell(ii,i)
   if(plug%end_type.eq.'begin ') then
    call insert(temp,t_new,'after ')
   else
    call insert(temp,t_new,'before')
   end if
   temp=>t_new
  end do

!   call check_list_c(5,'up    ')

! The following plug the curve into the current and the other 
! node cells.
  pl_new=plug; pl_new%cv_nmb=new_curve
  pl_new%side_front='left  '; pl_new%side_behind='right '
  pl_new%end_type='begin '
  pl_new%outgoing_wave_number=pp_cell(ii,1)%wv_nb
  allocate(tp_new); tp_new%plug=pl_new

!   plugg=tempp%plug; adsss=tempp%address
  select case(plug%end_type)
   case('begin ')
    call insert(tp_insert,tp_new,'before')
    if(insert_order.eq.'clock   ') tp_insert=>tp_new
   case('end   ')
    call insert(other_tp_insert,other_tp_new,'after ')
    if(insert_order.eq.'counter ') other_tp_insert=>other_tp_new
  end select
  other_new_plug=other_plug; other_new_plug%cv_nmb=new_curve
  other_new_plug%side_front='right '
  other_new_plug%side_behind='left  '
  other_new_plug%end_type='end   '
  allocate(other_tp_new); other_tp_new%plug=other_new_plug
  call insert(other_tempp,other_tp_new,'after ')
   
!   call check_ring('counter ')
!   call check_list_c(5,'down  ')

! Update the solution in smooth, update the grid-map and preserve 
! conservation.
  t_new=>cvv(new_curve)%begin%next
  head_mark=t_new%address%idx; head_switch=.true.
  do while(head_switch)
   call update_near_node(t_new,tp_new)

!    call check_list_c(5,'down  ')
!    call check_ring('counter ')

!   call update_near_node(t_new,other_tp_new)
!    call check_list_c(5,'down  ')

   t_new=>t_new%next
   head_switch=(associated(t_new))
  end do

 end if
end do
 
! call check_list_c(6,'down  ')

select case(plug%end_type)
 case('begin '); tempp=>tempp%previous
 case('end   '); tempp=>tp_new
end select

! call check_list_c(5,'down  ')

! Delete the old discontinuity curve and its curve-plug.
call deletee(cvv(plug%cv_nmb))

! call check_list_c(5,'down  ')

end subroutine yawn_update


end module tools_4_connecting_n_triple