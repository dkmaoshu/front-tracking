module node_production_n_update_2
! This module is for production of node cell information for the
! new node cell and build the discontinuity curve connecting the 
! old and new node cells.

use variables_4_node_treat
! 'bs_var.f90'

public  connect_n_produce
private connect_old_n_new, produce_node_info 


contains


subroutine connect_n_produce
! Implement the production and connection.

implicit none

if(current_node_empty.eq.'no ') call connect_old_n_new

! call check_ring(1,'counter ')

call produce_node_info

end subroutine connect_n_produce


subroutine connect_old_n_new
! This subroutine creats a discontinuity curve connecting the
! old and newly produced node cells.

implicit none
type(curve_plug), pointer :: temppo, temppn
type(cv_plug_info) :: plug_o, plug_n !, plugg
type(state), dimension(-1:1,-1:1) :: u_f, u_b
integer :: cv_nmb, edge_o, edge_n
character*6, dimension(2) :: direction
! To tell on what direction the new node cell is viewed from 
! the old node cell.
integer :: i_dif, j_dif, i, j

! The following produces the discontinuity curve connecting the
! old and new node cells.
do cv_nmb=1,cvn
 if(cvv(cv_nmb)%status.eq.'asleep') then
  cvv(cv_nmb)%status='yawn  '
  cvv(cv_nmb)%cv_type='dblink  '
  cvv(cv_nmb)%begin_end=nd_nmb; cvv(cv_nmb)%end_end=nd_nmb_n
  cvv(cv_nmb)%total=0; cvv(cv_nmb)%wave=error_data
  call creat_cv(cv_nmb)
  exit
 end if
end do

! Determine on what direcion the new node cell is viewed from 
! the old node cell.
direction='      '
i_dif=x_nd_new-x_nd; j_dif=y_nd_new-y_nd
select case(i_dif)
 case(-1); direction(1)='west  '
 case(0);  direction(1)='none  '
 case(1);  direction(1)='east  '
end select
select case(j_dif)
 case(-1); direction(2)='south '
 case(0);  direction(2)='none  '
 case(1);  direction(2)='north '
end select

select case(direction(1))
 case('east  ')
  edge_n=4
  select case(direction(2))
   case('north '); edge_o=3
   case('none  '); edge_o=2
   case('south '); edge_o=1
  end select
 case('none  ')
  select case(direction(2))
   case('north '); edge_o=3; edge_n=1
   case('none  '); print*,'Something wrong!!!'; pause
   case('south '); edge_o=1; edge_n=3
  end select
 case('west  ')
  edge_n=2
  select case(direction(2))
   case('north '); edge_o=3
   case('none  '); edge_o=4
   case('south '); edge_o=1
  end select
end select 

! The following produces the plugs that plug the discontinuity
! curve on the old and new node cell.
plug_o%cv_nmb=cv_nmb
plug_o%end_type='begin '
plug_o%edge=edge_o
plug_o%side_front='left  '; plug_o%side_behind='right '

do i=-1,1; do j=-1, 1
 u_f(i,j)=tempp_n%plug%u_behind(i,j)
 u_b(i,j)=tempp%plug%u_front(i,j)
end do; end do
allocate(temppo)

do i=-1,1; do j=-1, 1
 plug_o%u_front(i,j)=u_f(i,j)
 plug_o%u_behind(i,j)=u_b(i,j)
end do; end do
temppo%plug=plug_o

! call check_ring(1,'clock   ')

! plugg=tempn%plug

plug_n%cv_nmb=cv_nmb
plug_n%end_type='end   '
plug_n%edge=edge_n
plug_n%side_front='right '; plug_n%side_behind='left  '
do i=-1,1; do j=-1, 1
 u_f(i,j)=tempn%plug%u_behind(i,j)
 u_b(i,j)=tempn_n%plug%u_front(i,j)
end do; end do
allocate(temppn)
do i=-1,1; do j=-1, 1
 plug_n%u_front(i,j)=u_f(i,j)
 plug_n%u_behind(i,j)=u_b(i,j)
end do; end do
temppn%plug=plug_n

! Insert the curve-plugs corresponding to the newly produced 
! discontinuity curve into the olde and new node cells.

! call check_ring(1,'clock   ')

call insert(tempp,temppo,'after ')

! call check_ring(1,'counter ')
if(new_nodes.eq.1) head_mark=temppo%address%idx

call insert(tempn,temppn,'before ')

call renumber_plugs(nd_nmb_n)

end subroutine connect_old_n_new


subroutine produce_node_info
! This subroutine produces the node cell information for the
! new node cell.

implicit none
type(node_info) :: n_cell
integer :: i

n_cell%x_idx=x_nd_new; n_cell%y_idx=y_nd_new
n_cell%x_posi=xw_posi; n_cell%y_posi=yw_posi
do i=1,4
 n_cell%flux(i)=error_data
end do
ndd(nd_nmb_n)%n_cell=n_cell

! call check_ring(5,'counter ')
! call check_list_c(5,'up    ')

call evaluate_node_state(nd_nmb_n,ndd(nd_nmb_n)%n_cell%or_state)
! The ordinary cell-average in the new node cell is computed based
! on the configuration of the discontinuity curves meeting at the 
! node.

sum=sum-ndd(nd_nmb_n)%n_cell%or_state
ndd(nd_nmb)%n_cell%or_state=ndd(nd_nmb)%n_cell%or_state+sum
! The ordinary cell-average in the old node cell is updated based
! on preserving conservation.

call clean(ggd_cell(x_nd_new,y_nd_new))
ggd_cell(x_nd_new,y_nd_new)%region='node'
! Update the grid map.

end subroutine produce_node_info


end module node_production_n_update_2