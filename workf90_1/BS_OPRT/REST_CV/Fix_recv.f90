module fix_reversed_curves

use solu_comput
! 'solu_com.f90'

use node_cells
! 'nd_cells.f90'

implicit none


contains


subroutine fix_reversed_cv

implicit none
integer :: i

do i=1,cvn
 if(acvv(i)%status.eq.'awake ') then
  call fix_revers(i)
 end if
end do

! call check_list_ac(2,'up    ')

end subroutine fix_reversed_cv


subroutine fix_revers(curve_number)
! This subroutine fixes reversed curves. These curves consists of 
! a single 'xy'-type critical cell and this critical cell moves
! and thus its two ends reversed.

integer, intent(in) :: curve_number

type(geo_info) :: g_cell
type(node_info) :: n_cell_1, n_cell_2
type(curve_plug), pointer :: tempp
type(cv_plug_info) :: plug
integer :: i0, j0, i1, j1, i2, j2

if(acvv(curve_number)%total.ne.1) return
call clean_up_g_cell(g_cell)
call clean_up_n_cell(n_cell_1); call clean_up_n_cell(n_cell_2)
call clean_up_plug(plug)
if(acvv(curve_number)%cv_type.ne.'dblink') call error_message

! call delete(acvv(8))

g_cell=acvv(curve_number)%begin%next%g_cell
n_cell_1=ndd(acvv(curve_number)%begin_end)%n_cell
n_cell_2=ndd(acvv(curve_number)%end_end)%n_cell
i0=g_cell%x_idx; j0=g_cell%y_idx;
call find_neighboring_cell(i0,j0,g_cell%edge(1),i1,j1)
call find_neighboring_cell(i0,j0,g_cell%edge(2),i2,j2)
if(i1.ne.n_cell_1%x_idx.or.j1.ne.n_cell_1%y_idx) then
 if(i2.ne.n_cell_2%x_idx.or.j2.ne.n_cell_2%y_idx) then
  if(acvv(curve_number)%begin_end.ne.0) then
! Update the node-plug corresponding to the begin-end.
   call visit(acvv(curve_number)%begin_end,curve_number, &
              'begin ',tempp)
   plug=tempp%plug
   select case(plug%end_type)
    case('begin '); plug%end_type='end   '
	case('end   '); plug%end_type='begin '
   end select
   select case(plug%side_front)
    case('left  ')
	 plug%side_front='right '; plug%side_behind='left  '
	case('right ')
     plug%side_front='left  '; plug%side_behind='right '
    case default; call error_message
   end select
   plug%edge=neighbor_edge(g_cell%edge(2))
   tempp%plug=plug 
  end if
  if(acvv(curve_number)%end_end.ne.0) then
   call visit(acvv(curve_number)%end_end,curve_number, &
              'end   ',tempp)
   plug=tempp%plug
   select case(plug%end_type)
    case('begin '); plug%end_type='end   '
	case('end   '); plug%end_type='begin '
   end select
   select case(plug%side_front)
    case('left  ')
	 plug%side_front='right '; plug%side_behind='left  '
	case('right ')
     plug%side_front='left  '; plug%side_behind='right '
    case default; call error_message
   end select
   plug%edge=neighbor_edge(g_cell%edge(1))
   tempp%plug=plug
  end if
  call revers_auxi_curve(curve_number,'no ','yes','no ')
 end if
end if

end subroutine fix_revers


end module fix_reversed_curves