module produce_ghost_heads
! This module produces the ghost sections of curves at the two heads of a discontinuity
! curves. When the curve is 'awake ', the two ghost heads are just the 'begin%next' and
! 'end%previous' critical cells, and when the curve is 'yawn', the two ghost critical 
! cells are produced in someway.

use solution
! solution.f90'

implicit none


public  ghost_production, ghost_delete


contains


subroutine ghost_production
! This subroutine produces the ghost critical cells for all the discontinuity curves.

implicit none

integer :: i

do i=1, curves_number
 nullify(cvv(i)%begin_ghost); nullify(cvv(i)%end_ghost)
end do

do i=1, curves_number
 select case(cvv(i)%status)
  case('awake '); call gh_prod_awake(i)
  case('yawn  '); call gh_prod_yawn(i)
  case('asleep')
  case default; call error_message
 end select
end do


contains


subroutine gh_prod_awake(cv_nmb)

implicit none

integer, intent(in) :: cv_nmb

allocate(cvv(cv_nmb)%begin_ghost)
allocate(cvv(cv_nmb)%end_ghost)
cvv(cv_nmb)%begin_ghost=cvv(cv_nmb)%begin%next
cvv(cv_nmb)%end_ghost=cvv(cv_nmb)%eend%previous
cvv(cv_nmb)%message='orginal_awake  '

end subroutine gh_prod_awake


subroutine gh_prod_yawn(cv_nmb)

implicit none

integer, intent(in) :: cv_nmb

type(curve_plug), pointer :: tempp
integer :: head_mark, nd_nmb
logical :: head_switch
type(geo_info) :: g_cell
type(phy_info) :: p_cell

! Produce the 'begin_ghost' first.
! First find the adjacent node cell.
nd_nmb=cvv(cv_nmb)%begin_end

! Second find the corresponding curve-plug
tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)
 if(tempp%plug%cv_nmb.eq.cv_nmb.and.tempp%plug%end_type.eq.'begin ') exit 
 tempp=>tempp%next  
 head_switch=(tempp%address%idx.ne.head_mark)
end do
if(.not.head_switch) then
 print*, 'There must be something wrong.'; pause
end if

! Third compute the ghost 'g_cell'.
g_cell%dis_posi=0.0d0
select case(tempp%plug%edge)
 case(1)
  g_cell%x_idx=ndd(nd_nmb)%n_cell%x_idx; g_cell%y_idx=ndd(nd_nmb)%n_cell%y_idx-1
  g_cell%edge(1)=3; g_cell%edge(2)=1
 case(2)
  g_cell%x_idx=ndd(nd_nmb)%n_cell%x_idx+1; g_cell%y_idx=ndd(nd_nmb)%n_cell%y_idx
  g_cell%edge(1)=4; g_cell%edge(2)=2
 case(3)
  g_cell%x_idx=ndd(nd_nmb)%n_cell%x_idx; g_cell%y_idx=ndd(nd_nmb)%n_cell%y_idx+1
  g_cell%edge(1)=1; g_cell%edge(2)=3
 case(4)
  g_cell%x_idx=ndd(nd_nmb)%n_cell%x_idx-1; g_cell%y_idx=ndd(nd_nmb)%n_cell%y_idx
  g_cell%edge(1)=2; g_cell%edge(2)=4
end select 
call find_type_from_edges(g_cell)
call find_points_from_edges(g_cell)

! Fourth coompute the ghost 'p_cell'.
select case(tempp%plug%edge)
 case(1); p_cell%l_state=tempp%plug%u_front(0,-1); p_cell%r_state=tempp%plug%u_behind(0,-1)
 case(2); p_cell%l_state=tempp%plug%u_front(1,0);  p_cell%r_state=tempp%plug%u_behind(1,0)
 case(3); p_cell%l_state=tempp%plug%u_front(0,1);  p_cell%r_state=tempp%plug%u_behind(0,1)
 case(4); p_cell%l_state=tempp%plug%u_front(-1,0); p_cell%r_state=tempp%plug%u_behind(-1,0)
end select
p_cell%or_state=0.5d0*(p_cell%l_state+p_cell%r_state)

! Finally produce the 'ghost_begin' critical cell.
cvv(cv_nmb)%begin_ghost%g_cell=g_cell; cvv(cv_nmb)%begin_ghost%p_cell=p_cell


! Then produce the 'end_ghost'.
! First find the adjacent node cell.
nd_nmb=cvv(cv_nmb)%end_end

! Second find the corresponding curve-plug
tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)
 if(tempp%plug%cv_nmb.eq.cv_nmb.and.tempp%plug%end_type.eq.'end   ') exit 
 tempp=>tempp%next  
 head_switch=(tempp%address%idx.ne.head_mark)
end do
if(.not.head_switch) then
 print*, 'There must be something wrong.'; pause
end if

! Third compute the ghost 'g_cell'
g_cell%dis_posi=0.0d0
select case(tempp%plug%edge)
 case(1)
  g_cell%x_idx=ndd(nd_nmb)%n_cell%x_idx; g_cell%y_idx=ndd(nd_nmb)%n_cell%y_idx-1
  g_cell%edge(1)=1; g_cell%edge(2)=3
 case(2)
  g_cell%x_idx=ndd(nd_nmb)%n_cell%x_idx+1; g_cell%y_idx=ndd(nd_nmb)%n_cell%y_idx
  g_cell%edge(1)=2; g_cell%edge(2)=4
 case(3)
  g_cell%x_idx=ndd(nd_nmb)%n_cell%x_idx; g_cell%y_idx=ndd(nd_nmb)%n_cell%y_idx+1
  g_cell%edge(1)=3; g_cell%edge(2)=1 
 case(4)
  g_cell%x_idx=ndd(nd_nmb)%n_cell%x_idx-1; g_cell%y_idx=ndd(nd_nmb)%n_cell%y_idx
  g_cell%edge(1)=4; g_cell%edge(2)=2
end select 
call find_type_from_edges(g_cell)
call find_points_from_edges(g_cell)

! Finally coompute the ghost 'p_cell'.
select case(tempp%plug%edge)
 case(1); p_cell%r_state=tempp%plug%u_front(0,-1); p_cell%l_state=tempp%plug%u_behind(0,-1)
 case(2); p_cell%r_state=tempp%plug%u_front(1,0);  p_cell%l_state=tempp%plug%u_behind(1,0)
 case(3); p_cell%r_state=tempp%plug%u_front(0,1);  p_cell%l_state=tempp%plug%u_behind(0,1)
 case(4); p_cell%r_state=tempp%plug%u_front(-1,0); p_cell%l_state=tempp%plug%u_behind(-1,0)
end select
p_cell%or_state=0.5d0*(p_cell%l_state+p_cell%r_state)

! Finally produce the 'ghost_begin' critical cell.
cvv(cv_nmb)%end_ghost%g_cell=g_cell; cvv(cv_nmb)%end_ghost%p_cell=p_cell

cvv(cv_nmb)%message='orginal_ghost  '

end subroutine gh_prod_yawn


end subroutine ghost_production


subroutine ghost_delete

implicit none

integer :: i

do i=1, curves_number
 if(associated(cvv(i)%begin_ghost)) deallocate(cvv(i)%begin_ghost)
 if(associated(cvv(i)%end_ghost)) deallocate(cvv(i)%end_ghost)
end do

! call check_list_c(7,'down  ')

end subroutine ghost_delete


end module produce_ghost_heads