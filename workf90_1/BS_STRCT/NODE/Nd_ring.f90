module node_plug_ring
! This module describes curve plug-in ring of a node.

use physical_state
! 'burgers.f90'

use grid, ndn=>nodes_number
! 'grid.f90'

implicit none

type cv_plug_info
 integer :: cv_nmb
! The number of the discontinuity curve that is connected to the 
! node. 
  character*6 :: end_type
! The possible value for 'end_type' are 'begin' and 'end'.
  integer :: edge
! The cell-edge of the node cell the discontinuity curve plugs in.
 character*8 :: inout
! Variable indicating whether the plug-in curve is incoming or
! outgoing. This information is used in dealing with triple
! node cells
 integer :: outgoing_wave_number
! The wave number of outgoing discontinuity curves. This will be
! used in distributing conservation difference in dealing with
! movement of triple node cell.
 character*6 :: side_front, side_behind
! Two sides of the plug-in discontinuity curve.
 type(state), dimension(-1:1,-1:1) :: u_front, u_behind
! Local smooth solutions between the plug-in curve and the next
! and previous curve on the plug-in ring
end type cv_plug_info

type adss_plug
 integer :: nd_nmb
! The number of the node cell the plug is pluged.
 integer :: idx
! The identification number of the plug in the plug-ring.
end type adss_plug

type curve_plug
 type(adss_plug) :: address
! The identification of the plugging curve.
 type(cv_plug_info) :: plug
 type(curve_plug), pointer :: previous, next
end type curve_plug  

type curve_ring
! All the curves plug in the node cell form a circular doubly
! linked list (ring), which goes in an counter-clockwise way.
 type(curve_plug), pointer :: begin, eend
! The beginning and end of the curve plug-ring.
 integer :: total
! The total number of curve-plugs that plug in the node cell.
 integer :: top_idx
! The highest identification number of plugs in the ring.
 character*6 :: curves_type
! Indicating whether the curves pluging in the node are regular
! or auxiliary. The possible values for this variable are 'regula'
! and 'auxila'.

end type curve_ring

interface operator(.eq.)
 module procedure identify, compare_true
end interface

!interface clean_up
! module procedure clean_up_p, clean_up_address
!end interface

interface assignment(=)
 module procedure assign
end interface

interface operator(.ne.)
 module procedure compare_false
end interface

private identify, compare_true, compare_false
public  clean_up_address_plug, creat_ring, clean_up_plug


contains


subroutine clean_up_plug(plug)
! Set default value for a curve_plug.

implicit none

type(cv_plug_info), intent(out) :: plug
integer :: i, j

! real(8), dimension(2) :: pt

plug%cv_nmb=-1000
plug%end_type='      '
plug%edge=-1000
plug%side_front='      '; plug%side_behind='      '
do i=-1,1; do j=-1,1
 plug%u_front(i,j)=error_data; plug%u_behind(i,j)=error_data
end do; end do

! call find_corner(pt,2)
 
end subroutine clean_up_plug


function identify(plug1,plug2)	result(c)

implicit none
type(cv_plug_info), intent(in) :: plug1, plug2
logical :: c

c=.false.
if(plug1%cv_nmb.eq.plug2%cv_nmb) then
 if(plug1%end_type.eq.plug2%end_type) then
  if(plug1%edge.eq.plug2%edge) c=.true.
 end if
end if

end function identify


subroutine creat_ring(cv_rg)

implicit none
type(curve_ring), intent(out) :: cv_rg

allocate(cv_rg%begin); nullify(cv_rg%begin%previous)
nullify(cv_rg%begin%next); call clean_up_plug(cv_rg%begin%plug)
allocate(cv_rg%eend); nullify(cv_rg%eend%next)
nullify(cv_rg%eend%previous); call clean_up_plug(cv_rg%eend%plug)
cv_rg%total=0; cv_rg%top_idx=0; cv_rg%curves_type='regula'

end subroutine creat_ring


subroutine clean_up_address_plug(address)

implicit none
type(adss_plug), intent(out) :: address

address%nd_nmb=error_index
address%idx=error_index

end subroutine clean_up_address_plug


function compare_true(address1,address2) 

implicit none
type(adss_plug), intent(in) :: address1, address2
logical :: compare_true

compare_true=.false.
if(address1%nd_nmb.eq.address2%nd_nmb) then
 if(address1%idx.eq.address2%idx) compare_true=.true.
end if

end function compare_true


function compare_false(address1,address2) 

implicit none
type(adss_plug), intent(in) :: address1, address2
logical :: compare_false

compare_false=.false.
if(address1%nd_nmb.ne.address2%nd_nmb) then
 compare_false=.true.
else
 if(address1%idx.ne.address2%idx) compare_false=.true.
end if

end function compare_false


subroutine assign(address1,address2)

implicit none
type(adss_plug), intent(out) :: address1
type(adss_plug), intent(in) :: address2

address1%nd_nmb=address2%nd_nmb
address1%idx=address2%idx

end subroutine assign


end module node_plug_ring