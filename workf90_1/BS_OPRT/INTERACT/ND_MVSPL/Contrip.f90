module connect_and_triple
! This module treat the triple point nodes and complete 
! connections between discontinuity curves and the old node that
! are borken.

use fix_broken_connection
! 'fix_cnnt.f90'

use check_n_remove
! 'chk_rmv.f90'

use treatment_of_triples
! 'trt_trpl.f90'

use remove_broken_corners
! 'rm_b_crn.f90'

implicit none

public  connect_n_triple
private	connecting, new_geo_cell_n_area, remove_node, &
        check_triple, triple_treat, node_position


contains


subroutine connect_n_triple
! This subroutine performs the connection and triple node 
! treatment.

implicit none

character*8 :: tri_multiple

! integer :: iii

!sum=ndd(nd_nmb)%n_cell%or_state

! call check_ring('counter ')
! iii=nd_nmb

call check_triple(tri_multiple)

! call check_list_c(6,'up    ')

select case(tri_multiple)
 case('remove  '); call remove_node
 case('triple  ') 
  call remove_on_all_curves 
  call triple_treat
 case('multiple')
  
!   call check_list_c(4,'up    ')

  call remove_on_all_curves

!   call check_list_c(4,'up    ')

  call node_position

!   call check_list_c(4,'up    ')

  call connecting
  case default; print*, 'Something is wrong!!!'; pause
end select

! call check_list_c(2,'down  ')

end subroutine connect_n_triple


subroutine node_position
! This subroutine computes the node position in the old node cell.

implicit none

type(curve_plug), pointer :: tempp
type(cv_plug_info) :: plug
type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
character*3 :: if_empty
integer :: head_mark
logical :: head_switch
real(8), dimension(2,2) :: pt1, pt2
real(8), dimension(2) :: pt, posi
character*10 :: error_message
real(8) :: counter

nullify(tempp); nullify(temp)
tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
if_empty='no '; counter=1.0d0;  posi=0.0d0
do while(head_switch)
! The node position is computed as the mean-value of intersection
! points of all neighbored discontinuity curves (their extensions)
! if the intersection points lie in the node cell. Thus, we first
! compute these intersection points.
 plug=tempp%plug
 call locate_head(plug,'first ',temp,if_empty)
 if(if_empty.eq.'yes') exit
 g_cell=temp%g_cell
 call pick_point(g_cell,1,pt); pt1(1,:)=pt
 call pick_point(g_cell,2,pt); pt1(2,:)=pt
 
 plug=tempp%next%plug
 if(if_empty.eq.'yes') exit
 g_cell=temp%g_cell
 call pick_point(g_cell,1,pt); pt2(1,:)=pt
 call pick_point(g_cell,2,pt); pt2(2,:)=pt

 call ints_of_lines(pt1,pt2,pt,error_message)
 if(dabs(pt(1)).gt.0.5d0.or.dabs(pt(2)).gt.0.5d0) exit
 posi=(counter-1.0d0)/counter*posi+1.d0/counter*pt

 tempp=>tempp%next
 head_switch=(tempp%address%idx.ne.head_mark)
end do

ndd(nd_nmb)%n_cell%x_posi=posi(1)
ndd(nd_nmb)%n_cell%y_posi=posi(2)

end subroutine node_position


end module connect_and_triple