module check_n_remove
! This module checks the old node cell to see whether it is a
! triple or multiple node cell, or it should be removed.

use tools_4_connecting_n_triple
! 'tools_ct.f90'

use tools_used_in_triple
! 'trip_bg.f90'

implicit none

public  check_triple, remove_node
private	between_finder, triple_multiple


contains


subroutine check_triple(tri_multiple)
! This subroutine checks if node cell under concern is a triple
! node.

implicit none

character*8, intent(out) :: tri_multiple

type(curve_plug), pointer :: tempp, tempp_yawn
type(cv_plug_info) :: plug
integer :: head_mark, i
logical :: head_switch
type(curve_map), dimension(wave_number) :: curves_map

! type(cv_plug_info) :: pllg
! logical :: xxx

tri_multiple='        '
if(new_nodes.ne.1.and.new_nodes.ne.2) then
 tri_multiple='multiple'; return
end if
! If more than two new nodes are produced in splitting,	the node 
! cell must be multiple.

! First, find the 'yawn' discontinuity curve connecting the old
! and new node cells.

! xxx=(associated(ndd(nd_nmb)%cv_rg%begin%next))

if(.not.associated(ndd(nd_nmb)%cv_rg%begin%next)) then
 tri_multiple='remove  '; return
end if
 
tempp=>ndd(nd_nmb)%cv_rg%begin%next
do while(associated(tempp))
 plug=tempp%plug
 if(cvv(plug%cv_nmb)%status.eq.'yawn  ') then
  if(cvv(plug%cv_nmb)%begin_end.eq.nd_nmb_n) exit 
  if(cvv(plug%cv_nmb)%end_end.eq.nd_nmb_n) exit
 end if
 tempp=>tempp%next
end do
tempp_yawn=>tempp

! pllg=tempp_yawn%plug

! Then form the 'curves_map', which contains the information
! of the rest of discontinuity curves plug-in in the old node
! cell and is necessary for triple-multiple check.
do i=1,wave_number
 curves_map(i)%facing='no      '
 curves_map(i)%status='asleep'
end do
i=0
head_mark=tempp%address%idx; head_switch=.true.
tempp=>tempp%next

plug=tempp%plug
if(cvv(plug%cv_nmb)%status.eq.'yawn  ') then
 if(ndd(nd_nmb)%cv_rg%total.eq.2) then
  tempp%plug%inout='incoming'; tri_multiple='triple  '
 else
  tri_multiple='multiple'
 end if
 return
end if
! If there is an another 'yawn' status curve pluging in	the old 
! node cell, then if the total plug-in curves are two, it is
! 'triple', otherwise it is multiple.

do while(head_switch)
 i=i+1; if(i.gt.wave_number) exit
 tempp%plug%inout='outgoing'; plug=tempp%plug
 curves_map(i)%status=cvv(plug%cv_nmb)%status
 select case(cvv(plug%cv_nmb)%status)
  case('awake ')
   select case(wave_facing(cvv(plug%cv_nmb)%wave)) 
    case('left  ')
     if(cvv(plug%cv_nmb)%begin_end.eq.nd_nmb) then
      curves_map(i)%facing='next    '
     else 
      curves_map(i)%facing='previous'
     end if
    case('right ')
     if(cvv(plug%cv_nmb)%begin_end.eq.nd_nmb) then
      curves_map(i)%facing='previous'
     else
      curves_map(i)%facing='next    '
     end if
    case('no    ')
     curves_map(i)%facing='no       '
    case default; print*, 'Something is wrong!!!'; pause
   end select
  case('yawn  ')
   tri_multiple='multiple'; return
  case default; print*, 'Something  is wrong!!!'; pause
 end select
 tempp=>tempp%next
 head_switch=(tempp%address%idx.ne.head_mark)
end do

if(i.gt.wave_number) then
 tri_multiple='multiple'; return
end if
! If the incoming discontinuity curves are more than the allowed
! wave number the node cell must be multiple.

call triple_multiple(curves_map,tri_multiple)
if(tri_multiple.eq.'triple  ') then
 tempp=>tempp_yawn%next; i=1 !; pllg=tempp%plug
 do while(tempp%plug%inout.eq.'outgoing') 
  tempp%plug%outgoing_wave_number=curves_map(i)%number
  i=i+1; tempp=>tempp%next
 end do
end if

end subroutine check_triple


subroutine remove_node
! It may happen that the node cell itself moves to a new grid
! cell, in which case all the curve-plugs moves to plug in the
! new node cell. This subroutine deals with this case.

implicit none

type(curve_plug), pointer :: tempp, tempp_n, tempp_nbs
type(node_info) ::no_cell, nn_cell
type(cv_plug_info) :: plug, plug_n
type(state) :: value_f, value_b
integer :: io, jo, in, jn, edge, i !, j
character*3 :: if_between
logical :: head_switch, distributing
integer :: head_mark, number_of_curves
type(cv_plug_info), dimension(5) :: curve_plugs
type(curve_map), dimension(wave_number) :: curves_map
character*8 :: tri_multiple
!type(state), dimension(:), allocatable :: dist_conv
!type(state), dimension(:,:), allocatable :: dist_cov
type(critical_cell), pointer :: temp

sum=ndd(nd_nmb)%n_cell%or_state

! First find the edge between the old and new node cells. The 
! edge is viewed from the new node.
no_cell=ndd(nd_nmb)%n_cell
io=no_cell%x_idx; jo=no_cell%y_idx
nn_cell=ndd(nd_nmb_n)%n_cell
in=nn_cell%x_idx; jn=nn_cell%y_idx
edge=error_index
if(iabs(in-io)+iabs(jn-jo).ne.1) then
 print*, 'Something is wrong!!!'; pause
end if
if(in.gt.io) edge=4
if(in.lt.io) edge=2
if(jn.gt.jo) edge=1
if(jn.lt.jo) edge=3

! Next, find the two curves between which the grid cell of the
! old node cell locates.
tempp=>ndd(nd_nmb_n)%cv_rg%begin%next; plug=tempp%plug
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)
 tempp_n=>tempp%next; plug_n=tempp_n%plug
 call between_finder(plug,plug_n,edge,if_between)
 if(if_between.eq.'yes') exit
 tempp=>tempp_n; plug=tempp%plug
 head_switch=(tempp%address%idx.ne.head_mark)
end do
if(.not.head_switch) then
 print*, 'Something is wrong!!!'; pause
end if

! Finally, delete the old node cell, update the numerical
! solution and grid-map, and preserve the conservation.
value_f=plug%u_front(io-in,jo-jn)
value_b=plug_n%u_behind(io-in,jo-jn)
! Delete the node cell.
call delete_node(ndd(nd_nmb))
! Update grid-map
call clean(ggd_cell(io,jo)); uu(io,jo)=value_f
! Preserve the conservation.
sum=sum-value_f

! call check_ring(2,'counter ')

! The conservation difference 'sum' should be distributed to the 
! rest of the discontinuity curves.

! First, the number of the rest discontinuity curves should be 
! counted.
number_of_curves=0
do i=1,5; call clean_up_plug(curve_plugs(i)); end do

tempp_nbs=>ndd(nd_nmb_n)%cv_rg%begin%next
head_mark=tempp_nbs%address%idx; head_switch=.true.
do while(head_switch)
 tempp_nbs=>tempp_nbs%next;
 number_of_curves=number_of_curves+1
 curve_plugs(number_of_curves)=tempp_nbs%plug
 head_switch=(tempp_nbs%address%idx.ne.head_mark)
end do
call clean_up_plug(curve_plugs(number_of_curves))
number_of_curves=number_of_curves-2

! Determine whether distribution operation should be carried out.
distributing=.true.
do i=1,wave_number
 curves_map(i)%facing='no      '
 curves_map(i)%status='asleep'
end do
if(number_of_curves.gt.wave_number) then
 distributing=.false.
else
 do	i=1, number_of_curves
  select case(cvv(curve_plugs(i)%cv_nmb)%status)
   case('yawn  ') 
    distributing=.false.; exit
   case('awake ')
    curves_map(i)%status=cvv(curve_plugs(i)%cv_nmb)%status
    select case(wave_facing(cvv(curve_plugs(i)%cv_nmb)%wave)) 
     case('left  ')
      if(cvv(curve_plugs(i)%cv_nmb)%begin_end.eq.nd_nmb_n) then
       curves_map(i)%facing='next    '
      else 
       curves_map(i)%facing='previous'
      end if
     case('right ')
      if(cvv(curve_plugs(i)%cv_nmb)%begin_end.eq.nd_nmb_n) then
       curves_map(i)%facing='previous'
      else
       curves_map(i)%facing='next    '
      end if
     case('no    ')
      curves_map(i)%facing='no       '
     case default; print*, 'Something is wrong!!!'; pause
    end select
   case default; call error_message
  end select
 end do
 call triple_multiple(curves_map,tri_multiple)
! select case(tri_multiple)
!  case('multiple'); distributing=.false.
!  case('triple	'); distributing=.true.
!  case default; call error_message
! end select
 
 if(tri_multiple.eq.'triple  ') then
  distributing=.true.
 else
  if(tri_multiple.eq.'mulitiple') then
   distributing=.false.
  else
   call error_message
  end if
 end if

end if
 

! !!!!!!! INTERIM !!!!!!!
! Finally, distribute the conservation difference.
if(.not.distributing) then
 sum=sum+ndd(nd_nmb_n)%n_cell%or_state
 ndd(nd_nmb_n)%n_cell%or_state=sum
else
 sum=sum/dble(number_of_curves)
 do i=1, number_of_curves
  select case(curve_plugs(i)%end_type)
   case('begin '); temp=>cvv(curve_plugs(i)%cv_nmb)%begin%next
   case('end   '); temp=>cvv(curve_plugs(i)%cv_nmb)%eend%previous
   case default; call error_message
  end select
  temp%p_cell%or_state=temp%p_cell%or_state+sum
 end do
end if
!!!!!!!!!!!!!!!!! 

call change_node_number(nd_nmb_n,nd_nmb)

end subroutine remove_node


subroutine between_finder(plug,plug_n,edge,if_between)
! This subroutine finds the two discontinuity curves between
! which the old node cell locates.

type(cv_plug_info), intent(in) :: plug, plug_n
integer, intent(in) :: edge
character*3, intent(out) :: if_between

integer :: edge_plug, edge_plug_n, edge_between

if_between='no '
edge_plug=plug%edge
edge_plug_n=plug_n%edge
if(edge_plug.eq.edge.or.edge_plug_n.eq.edge) call error_message

if(edge_plug.eq.edge_plug_n) then
 return
else
 edge_between=cycle_4(edge_plug+1)
 do while(edge_between.ne.edge_plug)
  if(edge_between.eq.edge_plug_n) then
   return
  else
   if(edge_between.eq.edge) then
    if_between='yes'; return
   else
    edge_between=cycle_4(edge_between+1)
   end if
  end if
 end do
end if

end subroutine between_finder


end module check_n_remove