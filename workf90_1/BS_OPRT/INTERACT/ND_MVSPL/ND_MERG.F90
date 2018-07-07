module node_mergence
! This module is for the mergence of two node cells. After the 
! splitting and moving of nodes, two nodes may locate in the same
! grid cell. In this case, the they have to merge into one node.
! This case does happen, although seldom.

use solution
! 'solution.f90'

implicit none


contains


subroutine node_merge

implicit none

type(node_info) :: n_cell, n_cell_l
integer :: i, j

do i=1, nodes_number-1
 if(ndd(i)%status.eq.'awake ') then
  n_cell=ndd(i)%n_cell
 else
  cycle
 end if
 do j=i+1, nodes_number
  if(ndd(j)%status.eq.'awake ') then
   n_cell_l=ndd(j)%n_cell
  else
   cycle
  end if
  if(n_cell%x_idx.eq.n_cell_l%x_idx.and.n_cell%y_idx.eq.n_cell_l%y_idx) then
   call merge_the_two_nodes(i,j)
   exit
  end if
 end do
end do          


contains


subroutine merge_the_two_nodes(first,second)
! This subroutine implement the mergence.

implicit none

integer, intent(in) :: first, second

type(node_info) :: n_cell
type(curve_plug), pointer :: tempp, tempp_n, tempp_x, current
type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
integer :: i, single, head_mark, head_mark_p
type(cv_plug_info) :: plug
type(state) :: difference
logical :: if_delete, head_switch, head_switch_p, if_beginning
character*6 :: side

! First compute the 'node_info' of the merged node cell.
n_cell%x_idx=ndd(first)%n_cell%x_idx
n_cell%y_idx=ndd(first)%n_cell%y_idx
n_cell%x_posi=error_data; n_cell%y_posi=error_data
do i=1, 4
 n_cell%flux(i)=error_data
end do

! The ordinary cell-averages of the merged node is taken as the
! sum of the ordinary cell-averages of the two to-be-merged 
! node cells.
n_cell%or_state=ndd(first)%n_cell%or_state+ndd(second)%n_cell%or_state

! Then the twice-counted cell-averages of the disappeared critical
! cell should be deducted; meanwhile, the 'yawn' double-pluging curve and the corresponding
! plugs in the the two node cells should be deleted.
tempp=>ndd(first)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.; if_beginning=.true.
do while(head_switch)

 plug=tempp%plug; if_delete=.false.
 if(cvv(plug%cv_nmb)%status.eq.'yawn  ') then

  select case(plug%end_type)
   case('begin '); if(cvv(plug%cv_nmb)%end_end.ne.second) call error_message
   case('end   '); if(cvv(plug%cv_nmb)%begin_end.ne.second) call error_message
   case default; call error_message
  end select
  
! First deduct the double-counted cell-averages of the disappeared critical cell.
  select case(cvv(plug%cv_nmb)%message)
   case('orginal_awake  ')
    select case(plug%end_type)
     case('begin ')
      n_cell%or_state=n_cell%or_state-cvv(plug%cv_nmb)%begin_ghost%p_cell%or_state
     case('end   ')
      n_cell%or_state=n_cell%or_state-cvv(plug%cv_nmb)%end_ghost%p_cell%or_state
     case default; call error_message
    end select
   case('orginal_ghost  ')
   case default; call error_message
  end select
! Delete the disappeared curve.
  call deletee(cvv(plug%cv_nmb))

! Find the curve-plug of the curve at the other end and then delete the two curve-plugs.
  select case(plug%end_type)
   case('begin '); call visit(second,plug%cv_nmb,'end   ',tempp_n) 
   case('end   '); call visit(second,plug%cv_nmb,'begin ',tempp_n)
   case default; call error_message
  end select
  
  call deletee(tempp,'before')
  call deletee(tempp_n,'after ') 

  if(if_beginning) head_mark=tempp%address%idx
 
!   call check_ring('counter ')
   
 end if

 tempp=>tempp%next  
 head_switch=(tempp%address%idx.ne.head_mark); if_beginning=.false.
end do

! The double-pluging curve should be deleted and the conservation
! difference should be compensated to the ordinary cell-averages
! of the new node cell.
tempp=>ndd(first)%cv_rg%begin%next
head_mark_p=tempp%address%idx; head_switch_p=.true.
do while(head_switch_p)

 plug=tempp%plug
 if(cvv(plug%cv_nmb)%status.eq.'awake ') then
  select case(plug%end_type)
   case('begin ')
    if(cvv(plug%cv_nmb)%end_end.eq.second) if_delete=.true.
   case('end   ')
    if(cvv(plug%cv_nmb)%begin_end.eq.second) if_delete=.true.
  end select
  
  if(if_delete) then
   nullify(temp)
   select case(plug%end_type) 
    case('begin '); temp=>cvv(plug%cv_nmb)%begin%next
    case('end   '); temp=>cvv(plug%cv_nmb)%eend%previous
    case default; call error_message
   end select
   head_mark=temp%address%idx; head_switch=.true.
    	     
   do while(associated(temp).and.head_switch) 
    g_cell=temp%g_cell; p_cell=temp%p_cell
    call find_single(g_cell,single); side=g_cell%point(single)
    select case(plug%end_type)
     case('begin '); call remove(temp,'after ',side,difference)
     case('end   '); call remove(temp,'before',side,difference)
    end select	 
    	  
    n_cell%or_state=n_cell%or_state+difference
     
    if(associated(temp).and.temp%g_cell%x_idx.ne.error_index) then 
     head_switch=(temp%address%idx.ne.head_mark)
    else
     exit
    end if
   end do

! Delete the alraedy empty discontinuity curve.
   call deletee(cvv(plug%cv_nmb))

! Find the curve-plug of the curve at the other end and then delete the two curve-plugs.
   select case(plug%end_type)
    case('begin '); call visit(second,plug%cv_nmb,'end   ',tempp_n) 
    case('end   '); call visit(second,plug%cv_nmb,'begin ',tempp_n)
    case default; call error_message
   end select
   call deletee(tempp,'before')
   call deletee(tempp_n,'after ') 
   exit

  end if
 end if
  
 tempp=>tempp%next  
 head_switch_p=(tempp%address%idx.ne.head_mark_p)
end do

! Finally, plug all the curve-plugs of the second node cell onto the first one and then
! delete the second node cell.

head_mark=tempp_n%address%idx; head_switch=.true.; current=>tempp
do while(head_switch)
 allocate(tempp_x)
 tempp_x=tempp_n
 tempp_x%address%nd_nmb=first; tempp_x%address%idx=tempp_x%address%idx+10
 call insert(current,tempp_x,'after ')
 current=>tempp_x
 tempp_n=>tempp_n%next
 head_switch=(tempp_n%address%idx.ne.head_mark)
end do

call deletee(ndd(second))

! call check_ring('counter ')

ndd(first)%n_cell=n_cell
! ndd(second)%status='asleep'
 

end subroutine merge_the_two_nodes


end subroutine node_merge


end module node_mergence