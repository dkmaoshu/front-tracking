module variables_4_node_treat
! This module contains the basic variables and corresponding tools 
! used in the treatment of node move and splitting.

use solution
! 'solution.f90'

use polygon_intersection
! 'polyin.f90'

implicit none

integer :: nd_nmb, nd_nmb_n, new_nodes
! The numbers of the node cell being treated and new node cells
! produced in splitting and the total number of the new node cells.
integer :: x_nd, y_nd, x_nd_new=-1000, y_nd_new=-1000
! The x- and y-indexes of the node cell under concern and the
! will-be node cell.
real(8) :: xw_posi, yw_posi
! The x- and y-coordinates of the node point in the will-be node
! cell.
type(state) :: sum
! The total amount of physical quantity(ies) involved in the
! production. This variable is used to control the conservation.
character*3 :: if_merge
! Parameter telling whether there is merge of discontinuity curves.
character*6 :: current_merged, next_merged
! Indicate the critical cells on the current and next curves that
! will be merged.
character*3 :: current_node_empty
! Tell whether all the curve-plugs of the current node have move
! out.
type(curve_plug), pointer :: tempp, tempp_n, tempn, tempn_n
! Pointers 'tempp' and 'tempp_n' point first to the two detected
! curve-plugs that will merge outside of the old node cell. After
! the production of the new node cell and the update of 
! curve-plugs, they point to the two curve-plugs that form  the 
! 'open-mouth' in the plug-ring of the old node cell, in which a 
! curve-plug of the discontinuity curve connecting the old and 
! new node cell will be placed.
! Pointers 'tempn' and 'tempn_n' point to the two curve-plugs that
! form the 'open-mouth' in the plug-ring of the new node cell.
integer :: head_mark
logical :: head_switch


public nd_nmb, nd_nmb_n, x_nd, y_nd, x_nd_new, y_nd_new, &
       xw_posi, yw_posi, current_merged, next_merged, &
	   tempp, tempp_n, head_mark, head_switch, locate_head, &
	   find_node, find_edgec, evaluate_node_state, sum, &
       find_plug_posi, find_node_number, update_near_node


contains


subroutine find_node(g_cell,gn_cell,dif1,dif2,xw_posi,yw_posi)
! This subroutine finds the will-be node point coordinates.

implicit none
type(geo_info), intent(in) :: g_cell, gn_cell
real(8), dimension(2) :: dif1, dif2
real(8), intent(out) :: xw_posi, yw_posi

real(8), dimension(2) :: pt1, pt2, pt3, pt4, pt
real(8), dimension(2,2) :: pts1, pts2
character*10 :: err_msg

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)
pts1(1,:)=pt1+dif1; pts1(2,:)=pt2+dif1
call pick_point(gn_cell,1,pt3)
call pick_point(gn_cell,2,pt4)
pts2(1,:)=pt3+dif2; pts2(2,:)=pt4+dif2
call ints_of_lines(pts1,pts2,pt,err_msg)

select case(err_msg)
 case('ok        ')
  xw_posi=pt(1); yw_posi=pt(2)
 case('short     ')
  xw_posi=error_data; yw_posi=error_data
  return
 case('parallel  ')
  call error_message; pause
end select

end subroutine find_node


subroutine find_edgec(g_cell,plug,ii,jj,edge_c)
! When the head critical cell with 'g_cell' will be moved, it must 
! be located in a diagonal grid cell and there is a grid cell 
! between the diagonal grid cell and the node cell. This subroutine 
! finds the edge of the node cell adjacent to the grid cell in 
! between.

implicit none
type(geo_info), intent(in) :: g_cell
! The head critical cell.
type(cv_plug_info) :: plug
! The curve plug information.
integer, intent(out) :: edge_c
! The common edge to be found, which is viewed in the node cell.
integer, intent(out) :: ii, jj
! The coordinates of the grid cell in between.

select case(plug%edge)
 case(1,3)
  ii=g_cell%x_idx; jj=y_nd   
 case(2,4)
  ii=x_nd; jj=g_cell%y_idx
end select

select case(ii-x_nd)
 case(1); edge_c=2
 case(-1); edge_c=4
 case(0)
  select case(jj-y_nd)
   case(1); edge_c=3
   case(-1); edge_c=1
  end select
end select 

end subroutine find_edgec


subroutine find_plug_posi(g_cell,edge_c,posi,error_message)
! This subroutine finds the could-be plug-in discontinuity position
! for a head critical cell cell with broken connection with the
! node cell.

implicit none
type(geo_info), intent(in) :: g_cell
integer, intent(in) :: edge_c
real(8), intent(out) :: posi
character*10, intent(out) :: error_message

real(8), dimension(2) :: pt1, pt2, pt3, pt4, pt
real(8), dimension(2,2) :: pts1, pts2
real(8), dimension(2) :: dif
integer :: ii, jj, edge_n

error_message='          '
select case(edge_c)
 case(1); ii=x_nd; jj=y_nd-1
 case(2); ii=x_nd+1; jj=y_nd
 case(3); ii=x_nd; jj=y_nd+1
 case(4); ii=x_nd-1; jj=y_nd
end select
dif(1)=g_cell%x_idx-ii; dif(2)=g_cell%y_idx-jj

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)
pt1=pt1+dif; pt2=pt2+dif
pts1(1,:)=pt1; pts1(2,:)=pt2

edge_n=neighbor_edge(edge_c)
call pick_edge(edge_n,pt3,pt4)
pts2(1,:)=pt3; pts2(2,:)=pt4

call ints_of_lines(pts1,pts2,pt,error_message)

select case(edge_c)
 case(1,3); posi=pt(1)
 case(2,4); posi=pt(2)
end select

end subroutine find_plug_posi


subroutine locate_head(plug,position,temp,if_empty)
! This subroutine finds the head critical cell of pluging curve
! for given plug-in information.

implicit none
type(cv_plug_info), intent(in) :: plug
character*6, intent(in) :: position
type(critical_cell), pointer :: temp
character*3, intent(out) :: if_empty

nullify(temp); if_empty='no '

select case(plug%end_type)
 case('begin ')
  if(associated(cvv(plug%cv_nmb)%begin%next)) then
   select case(position)
    case('first ')
     temp=>cvv(plug%cv_nmb)%begin%next
    case('second')
     if(associated(cvv(plug%cv_nmb)%begin%next%next)) then
      temp=>cvv(plug%cv_nmb)%begin%next%next
     else
      return
     end if
    case('third ')
     if(associated(cvv(plug%cv_nmb)%begin%next%next)) then
      if(associated(cvv(plug%cv_nmb)%begin%next%next%next)) then 
       temp=>cvv(plug%cv_nmb)%begin%next%next%next
      else
       return
      end if
     end if
   end select
  else
   if_empty='yes'; return 	 	  	 	  	 	 
  end if
 case('end   ')
  if(associated(cvv(plug%cv_nmb)%eend%previous)) then
   select case(position)
    case('first ')
     temp=>cvv(plug%cv_nmb)%eend%previous
    case('second')
     if(associated(cvv(plug%cv_nmb)%eend%previous%previous)) then
      temp=>cvv(plug%cv_nmb)%eend%previous%previous
     else
      return
     end if
    case('third ')
     if(associated(cvv(plug%cv_nmb)%eend%previous%previous)) then
      if(associated(cvv(plug%cv_nmb)%eend%previous%previous% &
       previous)) then
  	   temp=>cvv(plug%cv_nmb)%eend%previous%previous%previous
      else
       return
      end if
     end if
   end select
  else
   if_empty='yes'; return 	 	  	 	  	 	 
  end if
end select 

end subroutine locate_head


subroutine find_node_number(nd_nmb)
! This program find node number for a new node cell.

implicit none
integer, intent(out) :: nd_nmb

integer :: i

nd_nmb=-1000
do i=1, ndn
 if(ndd(i)%status.eq.'asleep') then
  nd_nmb=i; exit
 end if
end do

end subroutine find_node_number


subroutine evaluate_node_state(nd_nmb,or_state)
! This subroutine computes the ordinary cell-average in a node
! cell based on the discontinuity curves meeting at the node.

implicit none

integer, intent(in) :: nd_nmb
! The number of the node cell whose ordinary cell-average is to
! be computed.
type(state), intent(out) :: or_state
! The computed ordinary cell-average.

type(curve_plug), pointer :: tempp
type(critical_cell), pointer :: temp
type(cv_plug_info), dimension(2) :: plug
! Curve-plugs of each two neighboring discontinuity curves in the
! plug-ring.
type(geo_info), dimension(2) :: g_cell
! The head critical cells of the two neighboring discontinuity
! curves.
type(node_info) :: n_cell, n_cell_n
! The current node cell and the node cell at the other end of
! a empty discontinuity curve.
integer :: head_mk
logical :: head_swt

type(polygon) :: p1, p2, p3
type(state), dimension(-1:1,-1:1) :: u_local
real(8) :: area, area_1, area_check !, aa(3)
integer :: i, j, l, ll, lll, x_dif, y_dif, first_edge, final_edge, medge 
character*3 :: if_empty, x_or_y
real(8), dimension(2) :: first_position, final_position, mposition
real(8), dimension(4,2) :: corners

or_state=0.d0; area_check=0.d0
n_cell=ndd(nd_nmb)%n_cell; call clean_up_n_cell(n_cell_n)
first_position=error_data; final_position=error_data; mposition=error_data
first_edge=error_index; final_edge=error_index; medge=error_index

! The second polygon is the new node cell.
p2%x(1)=-0.5d0; p2%y(1)=-0.5d0; p2%x(2)=-0.5d0; p2%y(2)=0.5d0
p2%x(3)=0.5d0; p2%y(3)=0.5d0; p2%x(4)=0.5d0; p2%y(4)=-0.5d0
p2%size=4

corners(1,1)=-0.5d0; corners(1,2)=-0.5d0
corners(2,1)=0.5d0;  corners(2,2)=-0.5d0
corners(3,1)=0.5d0;  corners(3,2)=0.5d0
corners(4,1)=-0.5d0; corners(4,2)=0.5d0

! The first polygon is the sector formed by each two neighbored
! discontinuity curves with its second vortex always being the
! node point.
p1%x=error_data; p1%y=error_data
p1%x(2)=n_cell%x_posi; p1%y(2)=n_cell%y_posi; p1%size=3
! ii=1

tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mk=tempp%address%idx; head_swt=.true.
do while(head_swt)

! The following forms the first polygon
 plug(1)=tempp%plug; plug(2)=tempp%next%plug; p1%size=3

! First download the local smooth solution.
 do i=-1, 1; do j=-1, 1
  u_local(i,j)=plug(1)%u_front(i,j)
 end do; end do

! Compute the first-and final-edges and the first-and final-positions.
do l=1, 2
 call locate_head(plug(l),'first ',temp,if_empty)
  select case(if_empty)

   case('no ')
! When the discontinuity curve is not empty, the discontinuity
! position of the head critical cell on the edges of the node
! cell is taken.
    g_cell(l)=temp%g_cell; medge=plug(l)%edge
    select case(plug(l)%edge)
     case(1); mposition(2)=-0.5d0; x_or_y='xxx'
     case(2); mposition(1)=0.5d0;  x_or_y='yyy'
     case(3); mposition(2)=0.5d0;  x_or_y='xxx'
     case(4); mposition(1)=-0.5d0; x_or_y='yyy'
    end select
    select case(plug(l)%end_type)
     case('begin ')
      select case(x_or_y)
       case('xxx'); mposition(1)=g_cell(l)%dis_posi(1)
       case('yyy'); mposition(2)=g_cell(l)%dis_posi(1)
      end select
     case('end   ')
      select case(x_or_y)
       case('xxx'); mposition(1)=g_cell(l)%dis_posi(2)
       case('yyy'); mposition(2)=g_cell(l)%dis_posi(2)
      end select
    end select

   case('yes')
! When the discontinuity curve is empty, the edge and the position are
! taken according to the position of the node cell at the other end
! of the discontinuity curve.
    if(plug(l)%end_type.eq.'begin ') then
     n_cell_n=ndd(cvv(plug(l)%cv_nmb)%end_end)%n_cell
    else 
     n_cell_n=ndd(cvv(plug(l)%cv_nmb)%begin_end)%n_cell
    end if
 
    x_dif=n_cell_n%x_idx-n_cell%x_idx
    y_dif=n_cell_n%y_idx-n_cell%y_idx
    select case(x_dif)
     case(-1); mposition(1)=-0.5d0
     case(0);  mposition(1)=0.0d0
     case(1);  mposition(1)=0.5d0
     case default; print*, 'Something wrong!!!'; pause
    end select
    select case(y_dif)
     case(-1); mposition(2)=-0.5d0
     case(0);  mposition(2)=0.0d0
     case(1);  mposition(2)=0.5d0
     case default; print*, 'Something wrong!!!'; pause
    end select
    if(iabs(x_dif)+iabs(y_dif).eq.0) call error_message
    if(iabs(x_dif)+iabs(y_dif).eq.1) then
    call find_edge_between(n_cell%x_idx,n_cell%y_idx,n_cell_n%x_idx, &
  	                       n_cell_n%y_idx,medge)
    else
     if(y_dif.gt.0) then
      medge=3; mposition(2)=0.5d0
     else
      medge=1; mposition(2)=-0.5d0
     end if
     if(x_dif.gt.0) then
      mposition(1)=0.4999999d0
     else
      mposition(1)=-0.4999999d0
     end if
    end if

  end select

  select case(l)
   case(1); first_edge=medge; first_position=mposition
   case(2); final_edge=medge; final_position=mposition
  end select

 end do

 ll=first_edge; p1%x(1)=first_position(1); p1%y(1)=first_position(2)
 area=0.d0
 do while(ll.ne.final_edge) 
  lll=cycle_4(ll+1)
  p1%x(3)=corners(lll,1); p1%y(3)=corners(lll,2)
  p3%x=error_data; p3%y=error_data; p3%size=0; area_1=0.d0
  call polyin(p1,p2,p3,'open  ',area_1) !; aa(ii)=area; ii=ii+1
  area=area+area_1
  p1%x(1)=p1%x(3); p1%y(1)=p1%y(3); ll=cycle_4(ll+1)
 end do
 p1%x(3)=final_position(1); p1%y(3)=final_position(2)
 p3%x=error_data; p3%y=error_data; p3%size=0; area_1=0.d0
 call polyin(p1,p2,p3,'open  ',area_1) !; aa(ii)=area; ii=ii+1
 area=area+area_1

 or_state=or_state+area*u_local(0,0)
 !write(*,*) ' Show something ', ue(st%f_region(l),dif_ix,dif_iy)
 area_check=area_check+area

 tempp=>tempp%next
 head_swt=(tempp%address%idx.ne.head_mk)

end do

if(dabs(area_check-1.0d0).gt.0.1d0) call error_message

! area=aa(1)+aa(2)+aa(3)

end subroutine evaluate_node_state


subroutine update_near_node(t_new,tp_new)
! This subroutine updates the numerical solution near a node cell
! when a new critical cell is inserted in a curve list near the 
! node either in completing connections of discontinuity curves 
! with the node or in dealing with triple node.

implicit none
type(critical_cell), pointer :: t_new
! The newly inserted critical cell (pointed by 't_new').
type(curve_plug), pointer :: tp_new
! The curve-plug of the discontinuity curve in which the critical
! cell is inserted.

type(curve_plug), pointer :: tp_nh
! The curve-plug of the discontinuity curve on which the possible 
! critical cell located in the same grid cell.
type(critical_cell), pointer :: temp_n
! The possible critical cell (pointed by 'temp_n') on the 
! neighboring discontinuity curve which is located in the same 
! grid cell.
type(cv_plug_info) :: pl_new, plug_n
type(geo_info) :: g_new, gn_cell
integer :: i0, j0
! The x- and y-indexes of the grid cell the critical cell under
! concern is located.
character*6, dimension(3) :: position
character*3 :: if_empty, neighbored
integer :: i
character*6 :: side_between, side_between_n, region, side_between_previous, side_between_next
character*6 :: side_sht
logical :: previous_neighbor_needed, next_neighbor_needed

previous_neighbor_needed=.false.; next_neighbor_needed=.false.
g_new=t_new%g_cell; i0=g_new%x_idx; j0=g_new%y_idx
call clean_up_g_cell(gn_cell)
position(1)='first '; position(2)='second'; position(3)='third '
region=ggd_cell(i0,j0)%region

! Find possible critical cell in the previous discontinuity curve 
! in the curve ring that is located in the same grid cell with
! 't_new'.
tp_nh=>tp_new%previous; nullify(temp_n); neighbored='no '
pl_new=tp_new%plug; plug_n=tp_nh%plug
select case(pl_new%end_type)
 case('begin '); side_between='right '
 case('end   '); side_between='left  '
 case default; print*, 'Something is wrong!!!'; pause
end select
select case(plug_n%end_type)
 case('begin '); side_between_n='left  '
 case('end   '); side_between_n='right '
 case default; print*, 'Something is wrong!!!'; pause
end select
do i=1,3
 call locate_head(plug_n,position(i),temp_n,if_empty)
 if(.not.associated(temp_n)) exit
 gn_cell=temp_n%g_cell
 if(gn_cell%x_idx.eq.i0.and.gn_cell%y_idx.eq.j0) then
  neighbored='yes'; exit
 end if
 call clean_up_g_cell(gn_cell)
end do
if(neighbored.eq.'yes') then
 select case(side_between)
  case('left  '); if(t_new%l_stk%cv_nmb.gt.0) call error_message
  case('right '); if(t_new%r_stk%cv_nmb.gt.0) call error_message
 end select
 call update_neighbor(t_new,temp_n,side_between,side_between_n)
 do i=1,4
  if(ggd_cell(i0,j0)%ccpt(i)%side.eq.side_between_n) then
   if(ggd_cell(i0,j0)%ccpt(i)%address.eq.temp_n%address) then
    side_sht=side_shift(side_between)
    if(g_new%point(i).eq.side_sht) then
     ggd_cell(i0,j0)%ccpt(i)%side=side_sht
!    select case(side_between)
!     case('left  '); ggd_cell(i0,j0)%ccpt(i)%side='right '
!     case('right '); ggd_cell(i0,j0)%ccpt(i)%side='left  '
!     case default; call error_message
!    end select
     ggd_cell(i0,j0)%ccpt(i)%address=t_new%address
    end if
   end if
  end if
 end do
else
 select case(region)
  
  case('smth  ','nd_smt')
   
   select case(side_between)
    case('left  ')
     if(t_new%l_stk%cv_nmb.eq.0) then
      call update_map(t_new,side_between)
     end if
    case('right ')
     if(t_new%r_stk%cv_nmb.eq.0) then
      call update_map(t_new,side_between)
     end if
   end select

  case('crit  ')
   
   side_between_previous=side_between
   previous_neighbor_needed=.true.
   
!   call insert_a_crit(side_between)

  case default; call error_message

 end select
end if

! call check_list_c(5,'down  ')

! Find possible critical cell in the next discontinuity curve
! in the curve ring that is located in the same grid cell with
! 't_new'.
tp_nh=>tp_new%next; nullify(temp_n); neighbored='no '
pl_new=tp_new%plug; plug_n=tp_nh%plug
select case(pl_new%end_type)
 case('begin '); side_between='left  '
 case('end   '); side_between='right '
 case default; print*, 'Something is wrong!!!'; pause
end select
select case(plug_n%end_type)
 case('begin '); side_between_n='right '
 case('end   '); side_between_n='left  '
 case default; print*, 'Something is wrong!!!'; pause
end select
do i=1,3
 call locate_head(plug_n,position(i),temp_n,if_empty)
 if(.not.associated(temp_n)) exit
 gn_cell=temp_n%g_cell
 if(gn_cell%x_idx.eq.i0.and.gn_cell%y_idx.eq.j0) then
  neighbored='yes'; exit
 end if
 call clean_up_g_cell(gn_cell)
end do
if(neighbored.eq.'yes') then
 select case(side_between)
  case('left  '); if(t_new%l_stk%cv_nmb.lt.0) call error_message
  case('right '); if(t_new%r_stk%cv_nmb.lt.0) call error_message
 end select
 call update_neighbor(t_new,temp_n,side_between,side_between_n)
 do i=1,4
  if(ggd_cell(i0,j0)%ccpt(i)%side.eq.side_between_n) then
   if(ggd_cell(i0,j0)%ccpt(i)%address.eq.temp_n%address) then
    side_sht=side_shift(side_between)
    if(g_new%point(i).eq.side_sht) then
     ggd_cell(i0,j0)%ccpt(i)%side=side_sht
!    select case(side_between)
!     case('left  '); ggd_cell(i0,j0)%ccpt(i)%side='right '
!     case('right '); ggd_cell(i0,j0)%ccpt(i)%side='left  '
!     case default; call error_message
!    end select
     ggd_cell(i0,j0)%ccpt(i)%address=t_new%address
    end if
   end if
  end if
 end do
else
 select case(region)
  
  case('smth  ','nd_smt')
   
   select case(side_between)
    case('left  ')
     if(t_new%l_stk%cv_nmb.eq.0) then
      call update_map(t_new,side_between)
     end if
    case('right ')
     if(t_new%r_stk%cv_nmb.eq.0) then
      call update_map(t_new,side_between)
     end if
   end select
   
  case('crit  ')

   side_between_next=side_between
   next_neighbor_needed=.true.
 !  call insert_a_crit(side_between)
    
   
  case default; call error_message 
   
 end select

end if

if(previous_neighbor_needed) call insert_a_crit(side_between_previous)
if(next_neighbor_needed) call insert_a_crit(side_between_next)


contains


subroutine insert_a_crit(side_between)

implicit none
character*6, intent(in) :: side_between

type(adss_info) :: address_stack
type(critical_cell), pointer :: temp_g, temp_gg
type(geo_info) :: g_cell_g, g_cell_gg
character*6 :: side, other_side, side_g
integer :: i

do i=1,4
 if(g_new%point(i).eq.side_between) then
  address_stack=ggd_cell(i0,j0)%ccpt(i)%address
  if(address_stack.ne.t_new%address) exit
 end if
end do
if(i.gt.4) return

call visit(address_stack,temp_gg)
g_cell_gg=temp_gg%g_cell
call which_side(g_new,g_cell_gg,side); side_g=side
call clean_up_g_cell(g_cell_g); nullify(temp_g)
do while(side_g.eq.side_between) 
 temp_g=>temp_gg; g_cell_g=g_cell_gg
 call which_side(g_cell_gg,g_new,other_side)
 select case(other_side)
  case('left  ')
   if(temp_g%l_stk%cv_nmb.gt.0) then
    address_stack=temp_g%l_stk
    call visit(address_stack,temp_gg)
	g_cell_gg=temp_gg%g_cell
	call which_side(g_new,g_cell_gg,side_g)
   else
    exit
   end if
  case('right ')
   if(temp_g%r_stk%cv_nmb.gt.0) then
    address_stack=temp_g%r_stk
    call visit(address_stack,temp_gg)
	g_cell_gg=temp_gg%g_cell
	call which_side(g_new,g_cell_gg,side_g)
   else
    exit
   end if
 end select
end do

if(side.eq.side_between) then
 select case(side)
  case('left  '); t_new%l_stk=temp_g%address
  case('right '); t_new%r_stk=temp_g%address
  case default; call error_message
 end select
 call which_side(g_cell_g,g_new,other_side)
 select case(other_side)
 case('left  '); temp_g%l_stk=t_new%address
 case('right '); temp_g%r_stk=t_new%address
 case default; call error_message
end select
else
 do i=1,4
  if(g_new%point(i).eq.side_between) then
   ggd_cell(i0,j0)%ccpt(i)%address=t_new%address
   ggd_cell(i0,j0)%ccpt(i)%side=side_between
  end if
 end do        
end if

end subroutine insert_a_crit


end subroutine update_near_node


end module variables_4_node_treat