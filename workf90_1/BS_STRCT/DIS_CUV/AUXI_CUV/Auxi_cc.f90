module auxiliary_critical_cells
! This module describes the data structure of auxlilary critiical
! cells.

use adss_info_cell
! 'adss_cell.f90'

use geo_info_cell
! 'gif_cell.f90'

use phy_info_cell
! 'pif_cell.f90'

use computational_information
! 'fif_cell.f90'

use critical_cells
! 'crit_cell.f90'

implicit none

type auxi_crit_cell

 type(adss_info) :: address
! The address of the auxiliary critical cell on the curve it 
! belongs, the 'address%idx' is counted from the 'begin'.

 type(adss_info) :: transfered
! The address of the critical cell to which the information in
! this auxiliary critical cell is transfered in the formation of
! discontinuity curves.

 type(geo_info) :: g_cell
! Geometrical information in auxiliary critical cell.

 type(phy_info) :: p_cell
! Physical information in auxiliary critical cell.

 type(flux_info) :: f_cell
! Flux information in auxiliary critical cell.

 type(comp_info) :: c_cell
! The information in the cell that is required in the computation,
! such as the physical states, discontinuity positions, and 
! numerical fluxes across the four edges on different levels in 
! Predictor-Corrector procedure.

 type(state), dimension(-1:1) :: ulx_exd, uly_exd, urx_exd,	ury_exd
! The data (or extrapolated data) of solution in the vicinity 
! of the auxiliary critical cell. These data will be used in 
! the reset of discontinuity curves when moves and splittings 
! occur.

 type(state) :: udiff_xl, udiff_yl, udiff_xr, udiff_yr
! Difference quotient in $x$- and $y$-directions on two sides.

 type(state) :: udiff_xl_f, udiff_yl_f, udiff_xr_f, udiff_yr_f
! Forward difference quotient in $x$- and $y$-directions on two sides.

 type(state) :: udiff_xl_b, udiff_yl_b, udiff_xr_b, udiff_yr_b
! Backward difference quotient in $x$- and $y$-directions on two sides.

real(8), dimension(3,2) :: stencil
! The stencil for interpolation used in reset procedure.
integer, dimension(3) :: stencil_disp
! To indicate the correspondence between the stencil points and discontinuity
! points.

 character*8 :: partner, partner_memory
! The partnership of neighbouring auxiliary critical cells.

 character*3 :: updated
! To tell whether the critical cell has been updated. This
! information is used in the 'reset' operation.

 character*3 :: corner
! Indicating if the critical cell is a corner of the curve.
 
 character*3 :: tangled
! Indicating if the critical cell is tangled. This information is
! used in 'reset'.

 type(auxi_crit_cell), pointer :: previous, next
! Neighboring information consists of two pointers pointing to 
! the previous and next critical cells of the critical cell on 
! the discontinuity curve
 
 type(adss_info) :: l_stk, r_stk
! Two pointers pointing to the left and right neighbors stacked 
! in the same grid.

 type(auxi_crit_cell), pointer :: p_nxt_dr, n_nxt_dr
! Two pointers pointing to its previous and next door neighbors.

 type(critical_cell), pointer :: temp

 type(state) :: memo

 real(8), dimension(2) :: memo_point

end type auxi_crit_cell


public  clean_up_auxi_crit


contains


subroutine clean_up_auxi_crit(tempa)

implicit none

type(auxi_crit_cell), pointer :: tempa

call clean_up_address(tempa%address)
call clean_up_g_cell(tempa%g_cell); call clean_up_p_cell(tempa%p_cell)
call clean_up_c_cell(tempa%c_cell); call clean_up_f_cell(tempa%f_cell)
call clean_up_address(tempa%l_stk); call clean_up_address(tempa%r_stk)
tempa%partner='single'; tempa%updated='no '

end subroutine clean_up_auxi_crit


subroutine find_edge_in_between_4_partnered(tempa,edge_between)

implicit none

type(auxi_crit_cell), pointer :: tempa
integer, intent(out) :: edge_between

type(geo_info) :: g_cell, g_partner

edge_between=error_index

select case(tempa%partner)
 case('previous'); g_partner=tempa%previous%g_cell
 case('next    '); g_partner=tempa%next%g_cell
 case('single  '); return
end select
g_cell=tempa%g_cell

if(g_partner%x_idx.gt.g_cell%x_idx) then
 edge_between=2
else
 if(g_partner%x_idx.lt.g_cell%x_idx) then
  edge_between=4
 else
  if(g_partner%y_idx.gt.g_cell%y_idx) then
   edge_between=3
  else
   if(g_partner%y_idx.lt.g_cell%y_idx) then
    edge_between=1
   else
    call error_message
   end if
  end if
 end if
end if

end subroutine find_edge_in_between_4_partnered


end module auxiliary_critical_cells