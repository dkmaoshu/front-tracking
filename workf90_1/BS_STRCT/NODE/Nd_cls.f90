module node_cells
! This module describes the data structure of the geometrical and 
! physical information in node cells.

use node_plug_ring
! 'nd_ring.f90'

use discontinuity_curves
! 'discontinuity_curve.f90'

type node_info
 integer :: x_idx, y_idx
! The x- and y-indexes of the node cell.
 real(8) :: x_posi, y_posi
 ! The x- and y-coordinates of the node point.
 type(state) :: or_state
! The orginary cell-average of the physical state in the node cell.
 type(state), dimension(4) :: flux
! The numerical fluxes on the four cell-edges.
end type node_info

type node_cell
 character*6 :: status
! Indicating if the node-cell is active.
 character*8 :: n_type
! Indicating the type of the node cell, whether it is 'triple' or
! 'multiple'
 type(node_info) :: n_cell
! The geometrical and physical information in the node cell.
 type(curve_ring) :: cv_rg
! The curve ring associated to the node cell.
end type node_cell

type(node_cell), dimension(ndn) :: nd
type(node_cell), dimension(:), allocatable :: ndd

interface visit
 module procedure visit_plug, visit_curve
end interface

interface insert
 module procedure insert_plug
end interface

interface deletee
 module procedure delete_plug, delete_node
end interface

private nd, visit_p, delete_plug, creat_ring, &
        insert_plug, visit_plug, visit_curve
public  ndd, give_nd, get_nd, clean_up_n_cell, visit, creat_node, &
        insert, curves_type_change, check_ring, deletee, &
		renumber_plugs, delete_node


contains


subroutine give_nd

implicit none
integer :: i

! real(8), dimension(2) :: pt

do i=1,ndn
 ndd(i)%status=nd(i)%status
 if(ndd(i)%status.eq.'awake') then
  ndd(i)%n_type=nd(i)%n_type
  ndd(i)%n_cell=nd(i)%n_cell
  ndd(i)%cv_rg%total=nd(i)%cv_rg%total
  ndd(i)%cv_rg%top_idx=nd(i)%cv_rg%top_idx
  ndd(i)%cv_rg%begin=>nd(i)%cv_rg%begin
  ndd(i)%cv_rg%eend=>nd(i)%cv_rg%eend
  ndd(i)%cv_rg%curves_type=nd(i)%cv_rg%curves_type
 end if
end do

! call find_corner(pt,2) 

end subroutine give_nd


subroutine get_nd

implicit none
integer :: i

do i=1,ndn
 nd(i)%status=ndd(i)%status
 if(ndd(i)%status.eq.'awake') then
  nd(i)%n_type=ndd(i)%n_type
  nd(i)%n_cell=ndd(i)%n_cell
  nd(i)%cv_rg%total=ndd(i)%cv_rg%total
  nd(i)%cv_rg%top_idx=ndd(i)%cv_rg%top_idx
  nd(i)%cv_rg%begin=>ndd(i)%cv_rg%begin
  nd(i)%cv_rg%eend=>ndd(i)%cv_rg%eend
  nd(i)%cv_rg%curves_type=ndd(i)%cv_rg%curves_type
 end if
end do 

end subroutine get_nd


subroutine check_ring(dir)
! This subroutine checks the curve-plug ring associated with a
! node cell.

implicit none

character*8, intent(in) :: dir
! The check direction, whether it's clockwise or counter-clockwise.
! The possible values for this variable are 'counter' and 'clock'.

integer :: nd_nmb
! The number of the node cell whose curve-plug ring is to be 
! checked.
type(curve_plug), pointer :: temp
type(cv_plug_info) :: plug_c, plug_p, plug_n
type(adss_plug) :: address_c, address_p, address_n
integer :: head_mark
logical :: head_switch

print*, 'The curve-plug ring of which node do you want to check?'
read(*,*) nd_nmb

if(ndd(nd_nmb)%status.eq.'asleep') then
 print*, 'The node is still sleeping.'; return
end if
if(.not.associated(ndd(nd_nmb)%cv_rg%begin)) then
 print*, 'The plug-ring is still empty yet.'; return
end if

nullify(temp)
if(dir.eq.'counter ') then
 if(associated(ndd(nd_nmb)%cv_rg%begin%next)) then
  temp=>ndd(nd_nmb)%cv_rg%begin%next
 else
  print*, 'The curve-plug ring is empty.'; return
 end if
else
 if(associated(ndd(nd_nmb)%cv_rg%eend%previous)) then
  temp=>ndd(nd_nmb)%cv_rg%eend%previous
 else
  print*, 'The curve-plug ring is empty.'; return
 end if
end if
head_mark=temp%address%idx; head_switch=.true.

do while(head_switch)
 call clean_up_plug(plug_c); call clean_up_plug(plug_p)
 call clean_up_plug(plug_n)

 plug_c=temp%plug; address_c=temp%address
 if(associated(temp%previous)) then
  plug_p=temp%previous%plug; address_p=temp%previous%address
 end if   
 if(associated(temp%next)) then
  plug_n=temp%next%plug; address_n=temp%next%address
 end if

 if(dir.eq.'counter ') then
  temp=>temp%next
 else
  temp=>temp%previous
 end if
 head_switch=(temp%address%idx.ne.head_mark)

end do

end subroutine check_ring


subroutine clean_up_n_cell(n_cell)
! Set default value for information in a node cell.

implicit none
type(node_info), intent(out) :: n_cell

n_cell%or_state%value=error_data
n_cell%x_idx=-1000; n_cell%y_idx=-1000

end subroutine clean_up_n_cell


subroutine visit_plug(nd_nmb,plug,tempp)
! This subroutine visits the curve-plug for given curve number
! 'cv_nmb', node number 'nd_nmb', and the end type of the plug.

implicit none

integer, intent(in) :: nd_nmb
type(cv_plug_info), intent(in) :: plug
type(curve_plug), pointer :: tempp

integer :: head_mark
logical :: head_switch

tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)
 if(tempp%plug.eq.plug) exit 
 tempp=>tempp%next  
 head_switch=(tempp%address%idx.ne.head_mark)
end do
if(.not.head_switch) then
 print*, 'There must be something wrong.'; pause
end if

end subroutine visit_plug


subroutine visit_curve(nd_nmb,cv_nmb,end_type,tempp)
! This subroutine visits the curve-plug for given curve number
! 'cv_nmb', node number 'nd_nmb', and the end type of the plug.

implicit none

integer, intent(in) :: nd_nmb, cv_nmb
character*6, intent(in) :: end_type
! Indicate whether the plug is a begin or an end.
type(curve_plug), pointer :: tempp

integer :: head_mark
logical :: head_switch

tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)
 if(tempp%plug%cv_nmb.eq.cv_nmb.and. &
    tempp%plug%end_type.eq.end_type) exit 
 tempp=>tempp%next  
 head_switch=(tempp%address%idx.ne.head_mark)
end do
if(.not.head_switch) then
 print*, 'There must be something wrong.'; pause
end if

end subroutine visit_curve


subroutine curves_type_change(curves_type)

implicit none
character*6, intent(in) :: curves_type

integer :: i

do i=1, ndn
 if(ndd(i)%status.eq.'awake ') then
  ndd(i)%cv_rg%curves_type=curves_type
 end if
end do

end subroutine curves_type_change


subroutine insert_plug(tempp,tempp_new,position)

implicit none
type(curve_plug), pointer :: tempp, tempp_new
character*6, intent(in) :: position

type(cv_plug_info) :: plug, plug_new
type(adss_plug) :: adss, adss_b, adss_e, adss_n
integer :: nd_nmb

adss=tempp%address; adss_n=adss; nd_nmb=adss%nd_nmb
ndd(nd_nmb)%cv_rg%total=ndd(nd_nmb)%cv_rg%total+1
ndd(nd_nmb)%cv_rg%top_idx=ndd(nd_nmb)%cv_rg%top_idx+1
adss_n%idx=ndd(nd_nmb)%cv_rg%top_idx
tempp_new%address=adss_n
plug=tempp%plug; plug_new=tempp_new%plug
nullify(tempp_new%previous); nullify(tempp_new%next)
call clean_up_address_plug(adss_b); call clean_up_address_plug(adss_e)

select case(position)
 case('before')
  if(associated(tempp%previous)) then
! The plug_ring is not empity.   
   tempp_new%previous=>tempp%previous
   tempp_new%next=>tempp
   tempp%previous%next=>tempp_new
   tempp%previous=>tempp_new
! The insert point may be the first plug in the ring; if so, 
! corresponding update should be done.
   adss_b=ndd(nd_nmb)%cv_rg%begin%next%address
   if(adss_b.eq.adss) ndd(nd_nmb)%cv_rg%begin%next=>tempp_new
  else
! The curve is still empty and 'tempp_new' is going to be the 
! first member of the curve. 'tempp' is the curve's 'begin'.
   tempp_new%previous=>tempp_new
   tempp_new%next=>tempp_new
   ndd(nd_nmb)%cv_rg%begin%next=>tempp_new
   ndd(nd_nmb)%cv_rg%eend%previous=>tempp_new
  end if

 case('after ')
  if(associated(tempp%next)) then
! The curve ring is not empty.  
   tempp_new%previous=>tempp
   tempp_new%next=>tempp%next
   tempp%next%previous=>tempp_new
   tempp%next=>tempp_new
! if the insert point may be the last plug in the ring; if so,
! corresponding update should be done.
   adss_e=ndd(nd_nmb)%cv_rg%eend%previous%address
   if(adss_e.eq.adss) ndd(nd_nmb)%cv_rg%eend%previous=>tempp_new  
  else
   tempp_new%previous=>tempp_new
   tempp_new%next=>tempp_new
   ndd(nd_nmb)%cv_rg%begin%next=>tempp_new
   ndd(nd_nmb)%cv_rg%eend%previous=>tempp_new
  end if
 case default; print*, 'Something wrong!!!'; pause
end select

end subroutine insert_plug


subroutine delete_plug(tempp,position)

implicit none
type(curve_plug), pointer :: tempp
! The pointer pointing to the will-be deleted plug.
character*6, intent(in) :: position
! Indicate whether the pointer will point to the plug previous
! or next to the deleted plug after the delection.

type(curve_plug), pointer :: current
type(adss_plug) :: adss, adss_p, adss_n, adss_b, adss_e
character*10 :: head_move

adss=tempp%address 
adss_p=tempp%previous%address; adss_n=tempp%next%address
if(position.eq.'before') then
 if(adss%idx.ne.adss_p%idx) then
  current=>tempp%previous
 else
  current=>ndd(adss%nd_nmb)%cv_rg%begin
 end if
else
 if(adss%idx.ne.adss_n%idx) then
  current=>tempp%next
 else
  current=>ndd(adss%nd_nmb)%cv_rg%eend
 end if
end if

tempp%previous%next=>tempp%next
tempp%next%previous=>tempp%previous

! The delete may delete the first or last plug, it may even 
! empty the ring. If these happen, corresponding updates should
! be done.

! Determine the cases.
adss_b=ndd(adss%nd_nmb)%cv_rg%begin%next%address
adss_e=ndd(adss%nd_nmb)%cv_rg%eend%previous%address
head_move='          '
if(adss_b%idx.ne.adss%idx) then
 if(adss_e%idx.ne.adss%idx) then
  head_move='no        '
 else
  head_move='eend_move '
 end if
else
 if(adss_e%idx.ne.adss%idx) then
  head_move='begin_move'
 else
  head_move='empty_ring'
 end if
end if
! Update the heads.
select case(head_move)
 case('no        ')
 case('begin_move')
  ndd(adss%nd_nmb)%cv_rg%begin%next=>tempp%next
 case('eend_move ')
  ndd(adss%nd_nmb)%cv_rg%eend%previous=>tempp%previous
 case('empty_ring')
  nullify(ndd(adss%nd_nmb)%cv_rg%begin%next)
  nullify(ndd(adss%nd_nmb)%cv_rg%eend%previous)
 case default; print*, 'Something is wrong'; pause
end select

deallocate(tempp)
if(head_move.ne.'empty_list') then
 tempp=>current
else
 nullify(tempp)
end if

ndd(adss%nd_nmb)%cv_rg%total=ndd(adss%nd_nmb)%cv_rg%total-1

end subroutine delete_plug


subroutine creat_node(nd_nmb)
! This subroutine creats a new node cell. The new node cell will
! be put in 'yawn' mode, a mode previous to 'awake' mode.

implicit none
integer, intent(out) :: nd_nmb

integer :: i
type(adss_plug) :: adss_b, adss_e

do i=1, ndn
 if(ndd(i)%status.eq.'asleep') then
  nd_nmb=i; exit
 end if
end do

ndd(nd_nmb)%status='yawn  '; ndd(nd_nmb)%n_type='        '
adss_b%nd_nmb=nd_nmb; adss_b%idx=-1
adss_e%nd_nmb=nd_nmb; adss_e%idx=-2

call creat_ring(ndd(nd_nmb)%cv_rg)
ndd(nd_nmb)%cv_rg%begin%address=adss_b
ndd(nd_nmb)%cv_rg%eend%address=adss_e

call clean_up_n_cell(ndd(nd_nmb)%n_cell)

end subroutine creat_node


subroutine waking_nodes
! This subrouitine wakes all the node cells in 'yawn  ' mode.

implicit none
integer :: i

do i=1, ndn
 if(ndd(i)%status.eq.'yawn  ') ndd(i)%status='awake '
end do

end subroutine waking_nodes


subroutine renumber_plugs(nd_nmb)
! This subroutine re-number plugs in the plug-ring of a node
! cell

implicit none
integer, intent(in) :: nd_nmb

type(curve_plug), pointer :: tempp
integer :: i, total

total=ndd(nd_nmb)%cv_rg%total; ndd(nd_nmb)%cv_rg%top_idx=total
if(total.gt.0) then
 tempp=>ndd(nd_nmb)%cv_rg%begin%next
 do i=1, total
  tempp%address%idx=i
  if(i.lt.total) tempp=>tempp%next
 end do 
end if

end subroutine renumber_plugs


subroutine delete_node(nddd)
! This subroutine deletes a node cell

implicit none
type(node_cell), intent(inout) :: nddd

type(curve_plug), pointer :: tempp

! Clean-up the plug-ring first
if(associated(nddd%cv_rg%begin%next)) then
 tempp=>nddd%cv_rg%begin%next
 do while(associated(tempp%next))
  call deletee(tempp,'after ')
!  call check_ring(1,'counter ')
 end do
end if
deallocate(nddd%cv_rg%begin)
deallocate(nddd%cv_rg%eend)
nddd%cv_rg%total=0
nddd%cv_rg%top_idx=0
nddd%cv_rg%curves_type='      '

! Clean-up the node cell information next.
call clean_up_n_cell(nddd%n_cell)

! Put the node cell into sleep finally.
nddd%status='asleep'
nddd%n_type='        '

end subroutine delete_node

 
end module node_cells