!This module is for the grid used for the computation. It contains the mesh size 'h', the
! the cell number 'nn', and the mesh ratio 'r'. 

module grid

implicit none

integer :: nxll, nyll, nxll_boundary, nyll_boundary
! The x- and y-indices of the low-left corner of the computational 
! region. 

integer :: nx, ny, nx_boundary, ny_boundary
! The numbers of grid cells in x- and y-directions.

real(8) :: x_width, y_width
! The width

real(8) :: h, r, rr=0.157300001d0
! Mesh size and mesh ratio.

real(8) :: current_time, final_time

integer :: total_step

real(8), parameter :: error_data=-1.0d4, x_length=1.0d0

integer, parameter :: error_index=-1000

integer, parameter :: curves_number=20
! The number of discontinuity curves used in the code.

integer, parameter :: nodes_number=10
! The number of node cells used in the code.

character*3 :: auxiliary_curve_type
! To indicate in which type the auxiliary discontinuity curves
! will be built.

character*12 :: order_of_curve_reset
! The order of the curve reset. If it is 'right_way', the first curve is reset first; if it
! is 'inverse_way', the last curve is reset first.

character*12 :: boundary_type_1, boundary_type_2
! The types of boundary. The first one indicates the particular test problem, and the second indicates the Thompson's nonreflecting treatment.
real(8) :: moved_time

public  set_grid, find_edge_between, find_neighboring_cell, &
        error_message, change_auxi_type, change_order_of_curve_reset, &
		moved_time


contains


subroutine set_grid

implicit none

h=dble(x_width/(nx-1))
! The mesh size is the width of grid cell edge. The cell is 
! assumed to be square, i.e. $h_x=h_y=h$.

end subroutine set_grid

subroutine change_auxi_type

implicit none

select case(auxiliary_curve_type)
 case('xx ')
  auxiliary_curve_type='yy '
 case('yy ')
  auxiliary_curve_type='xx '
 case default; print*, 'There must be something wrong!'; pause
end select

end subroutine change_auxi_type


subroutine change_order_of_curve_reset

implicit none

select case(order_of_curve_reset)
 case('right_way   '); order_of_curve_reset='inverse_way '
 case('inverse_way '); order_of_curve_reset='right_way   '
 case default; print*, 'Something must be wrong!'
end select

end subroutine change_order_of_curve_reset


subroutine error_message

print*, 'Something is wrong here!!!'
pause

end subroutine error_message


subroutine find_edge_between(i1,j1,i2,j2,edge_between)
! This subroutine finds the edge between two grid cells, the edge
! is viewed from the first grid cell

implicit none
integer, intent(in) :: i1, j1, i2, j2
! The x- and y-indexes of the two grid cells.
integer, intent(out) :: edge_between

edge_between=error_index
if(i1.ne.i2.and.j1.ne.j2) call error_message

if(i1.ne.i2) then
 if(i1.gt.i2) then
  edge_between=4
 else
  edge_between=2
 end if
else
 if(j1.ne.j2) then
  if(j1.gt.j2) then
   edge_between=1
  else
   edge_between=3
  end if
 else
  call error_message
 end if
end if

end subroutine find_edge_between


subroutine find_neighboring_cell(i0,j0,edge,i1,j1)
! Given the indexes of a grid cell and the across edge, find the
! indexes of the neighboring grid cell.

implicit none
integer, intent(in) :: i0, j0, edge
integer, intent(out) :: i1, j1

select case(edge)
 case(1); i1=i0; j1=j0+1
 case(2); i1=i0+1; j1=j0
 case(3); i1=i0; j1=j0+1
 case(4); i1=i0-1; j1=j0
 case default; call error_message
end select

end subroutine find_neighboring_cell


end module grid