module treatment_of_triples
! This module is to provide the treatment of triple node cells

use fix_triple_nodes
! 'fix_trip.f90'

implicit none

public triple_treat
private awake_curves, yawn_curve


contains


subroutine triple_treat
! This subroutine handles triple node-cells.

implicit none

type(cv_plug_info) :: plug
integer :: head_mark, cv_nmb
logical :: head_switch

sum=ndd(nd_nmb)%n_cell%or_state
ggd_cell(x_nd,y_nd)%region='nd_smt'

! call check_list_c(2,'down  ')
! call check_ring('counter ')

! Find the the 'yawn' status discontinuity curve connecting
! the old and new_node cells, which is actually the outgoing
! discontinuity curve.
tempp_outgoing=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp_outgoing%address%idx; head_switch=.true.
cv_nmb=-1000
do while(head_switch)
 plug=tempp_outgoing%plug; cv_nmb=plug%cv_nmb
 if(cvv(cv_nmb)%status.eq.'yawn  ') then
  if(cvv(cv_nmb)%begin_end.eq.nd_nmb_n)	exit
  if(cvv(cv_nmb)%end_end.eq.nd_nmb_n) exit
 end if
 cv_nmb=-1000
 tempp_outgoing=>tempp_outgoing%next
 head_switch=(tempp_outgoing%address%idx.ne.head_mark)
end do
if(cv_nmb.lt.0) then
 print*, 'Something is wrong!!!'; pause
end if

! The cases of 'awake'- and 'yawn'-status outgoing discontinuity
! curves are treated differentially.
cv_nmb=tempp_outgoing%next%plug%cv_nmb
select case(cvv(cv_nmb)%status)
 case('awake '); call awake_curves
 case('yawn  '); call yawn_curve
 case default; print*, 'Something is wron!!!'; pause
end select
 
! call check_list_c(2,'down  ')

end subroutine triple_treat


end module treatment_of_triples