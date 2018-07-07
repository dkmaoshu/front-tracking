module disjoints_fix
! This module is written for fixxing disjoints between head
! critical cells of discontinuity curves and the node cell. The 
! so-called disjoint is that the end discontinuity position does 
! not locate at the cell-edge between the two grid cells. What 
! the fix does is to recompute the problem discontinuity position 
! by extrapolation and reshape the critical cell based on the 
! recomputed position geometrically and physicallly. The 
! conservation difference caused then is compensated in the node 
! cell.

use solution
! 'solution.f90'

type(curve_plug), pointer :: tempp
type(critical_cell), pointer :: temp
type(cv_plug_info) :: plug
type(geo_info) :: g_cell, gn_cell
type(phy_info) :: p_cell
type(node_info) :: n_cell
type(state) :: difference
integer :: end_edge, edge_between, x_nd, y_nd, nd_nmb

public  fix_all_disjoints
private fix_disjoint, plug, g_cell, gn_cell, p_cell, difference, &
        temp, end_edge, edge_between, n_cell, fix_disjoints, &
		nd_nmb, x_nd, y_nd, tempp


contains 


subroutine fix_all_disjoints

implicit none

do nd_nmb=1,ndn
 if(ndd(nd_nmb)%status.eq.'awake ') call fix_disjoints
end do

end subroutine fix_all_disjoints


subroutine fix_disjoints
! This subroutine fixes disjoints around the node cell.

implicit none

logical :: head_switch
integer :: head_mark

n_cell=ndd(nd_nmb)%n_cell; x_nd=n_cell%x_idx; y_nd=n_cell%y_idx
tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)
 plug=tempp%plug; nullify(temp); edge_between=error_index
 if(cvv(plug%cv_nmb)%status.eq.'awake ') then
! Do only with 'awake ', i.e. non-empty curve.  
  if(plug%end_type.eq.'begin ') then
   if(associated(cvv(plug%cv_nmb)%begin%next)) then
    temp=>cvv(plug%cv_nmb)%begin%next
   end if
  else
   if(associated(cvv(plug%cv_nmb)%eend%previous)) then
    temp=>cvv(plug%cv_nmb)%eend%previous
   end if
  end if    
  if(associated(temp)) then 
   g_cell=temp%g_cell
   if(iabs(x_nd-g_cell%x_idx)+iabs(y_nd-g_cell%y_idx).eq.1) then
    call find_edge_between(g_cell%x_idx,g_cell%y_idx, &
                           n_cell%x_idx,n_cell%y_idx,edge_between)
    select case(plug%end_type)
     case('begin '); end_edge=g_cell%edge(1)
     case('end   '); end_edge=g_cell%edge(2)
	 case default; call error_message
    end select
    if(end_edge.ne.edge_between) call fix_disjoint
! There must be something wrong if a disjoint is found.
   end if
  end if
 end if

 tempp=>tempp%next
 head_switch=(tempp%address%idx.ne.head_mark)
end do 

end subroutine fix_disjoints


subroutine fix_disjoint
! This subroutine fixes disjoint between head_cells of 
! discontinuity curves and the node cell.

implicit none

type(node_info) :: n_cell
type(adss_info) :: address
real(8) :: posi
integer :: i0, j0, i

! integer :: ccv_nmb, ccv_mnb

n_cell=ndd(nd_nmb)%n_cell; difference=error_data
address=temp%address

select case(edge_between)
 case(1,3)
  select case(end_edge)
   case(2); posi=0.4999999d0
   case(4); posi=-0.4999999d0
   case default; call error_message
  end select
 case(2,4)
  select case(end_edge)
   case(1); posi=-0.4999999d0
   case(3); posi=0.4999999d0
  end select
end select   

select case(plug%end_type)
 case('begin ')
  g_cell%dis_posi(1)=posi; g_cell%edge(1)=edge_between
 case('end   ')
  g_cell%dis_posi(2)=posi; g_cell%edge(2)=edge_between
end select

call find_points_from_edges(g_cell)
call find_type_from_edges(g_cell)

temp%g_cell=g_cell
plug%edge=neighbor_edge(edge_between)
tempp%plug=plug

i0=g_cell%x_idx; j0=g_cell%y_idx

! ccv_nmb=temp%l_stk%cv_nmb
! ccv_mnb=temp%r_stk%cv_nmb

do i=1,4
 if(ggd_cell(i0,j0)%ccpt(i)%address.eq.address) then
  ggd_cell(i0,j0)%ccpt(i)%side=g_cell%point(i)
 end if
end do

!if(temp%l_stk%cv_nmb.le.0) then
! do i=1,4
!  if(g_cell%point(i).eq.'left  ') then
!   ggd_cell(i0,j0)%ccpt(i)%side='left  '
!  end if
! end do
!end if  
!if(temp%r_stk%cv_nmb.le.0) then
! do i=1,4
!  if(g_cell%point(i).eq.'right ') then
!   ggd_cell(i0,j0)%ccpt(i)%side='right '
!  end if
! end do
!end if  

end subroutine fix_disjoint
 

end module disjoints_fix
