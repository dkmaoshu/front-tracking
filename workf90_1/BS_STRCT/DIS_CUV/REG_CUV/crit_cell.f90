module critical_cells
! This module describes the data structure of critical cells.

use geo_info_cell
! 'gif_cell.f90'

use phy_info_cell
! 'phy_cell.f90'

use adss_info_cell
! 'adss_cell.f90'

implicit none

type critical_cell

 type(adss_info) :: address

 type(geo_info) :: g_cell
! Geometrical information in critical cell.

 type(phy_info) :: p_cell
! Physical information in critical cell.

 type(state) :: l_store, r_store
! For storing temporary information.

 character*6 :: tangled

 type(critical_cell), pointer :: previous, next
 type(adss_info) :: l_stk, r_stk
! Neighboring information consists of two pointers pointing to 
! the previous and next critical cells of the critical cell on 
! the discontinuity curve and two pointers pointing to the left
! and right neighbors stacked in the same grid.

 type(state) :: udiff_xl, udiff_yl, udiff_xr, udiff_yr
! Difference quotient in $x$- and $y$-directions on two sides.

 type(state) :: memo

 logical :: operation_memo
! To indicate whether an operation has been implemented on the critical 
! cell. 

end type critical_cell

public  clean_up_crit


contains


subroutine clean_up_crit(temp)

implicit none

type(critical_cell), pointer :: temp

! real(8), dimension(2) :: pt

call clean_up_address(temp%address)
call clean_up_g_cell(temp%g_cell); call clean_up_p_cell(temp%p_cell)
call clean_up_address(temp%l_stk); call clean_up_address(temp%r_stk)
 
! call find_corner(pt,2)

end subroutine clean_up_crit


subroutine update_neighbor(t_new,temp_n,side_between,side_between_n)
! This subroutine updates neighboring relation when a newly 
! produced critical cell is inserted in a grid cell which is 
! already occupied by an another critical cell.

implicit none
type(critical_cell), pointer :: t_new, temp_n
! The newly produced and already existed critical cells.
character*6, intent(in) :: side_between, side_between_n
! The side between the two critical cells viewed from the new
! and old critical cells.

type(phy_info) :: p_new, p_cell_n
type(state) :: state_between, dif

p_new=t_new%p_cell; p_cell_n=temp_n%p_cell; dif=0.d0
! Let corresponding pointers in each critical cell point to the
! other.
select case(side_between)
 case('left  ')
  t_new%l_stk=temp_n%address
  state_between=0.5d0*p_new%l_state
 case('right ')
  t_new%r_stk=temp_n%address
  state_between=0.5d0*p_new%r_state
 case default; print*, 'Something is wrong!!!'; pause
end select
select case(side_between_n)
 case('left  ')
  temp_n%l_stk=t_new%address
  state_between=state_between+0.5d0*p_cell_n%l_state
 case('right ')
  temp_n%r_stk=t_new%address
  state_between=state_between+0.5d0*p_cell_n%r_state
 case default; print*, 'Something is wrong!!!'; pause
end select

! Update the physcial state in between and the two ordinary
! states.
select case(side_between)
 case('left  ')
  p_new%l_state=state_between
 case('right ')
  p_new%r_state=state_between
end select
select case(side_between_n)
 case('left  ')
  dif=p_cell_n%l_state-state_between
  p_cell_n%l_state=state_between
 case('right ')   
  dif=p_cell_n%r_state-state_between
  p_cell_n%r_state=state_between
end select
p_cell_n%or_state=p_cell_n%or_state-dif
t_new%p_cell=p_new; temp_n%p_cell=p_cell_n
 
end subroutine update_neighbor


end module critical_cells