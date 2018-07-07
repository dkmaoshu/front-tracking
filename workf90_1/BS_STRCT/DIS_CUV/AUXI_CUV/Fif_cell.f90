module computational_information

use physical_state
! 'burgers.f90'

use grid
! 'grid.f90'

implicit none

type comp_info
 type(state), dimension(0:2) :: l_state, r_state, or_state
! Temporary storages for cell-averages on three time-levels in
! computation.
 type(state), dimension(4,2) :: flux
! Temporary storage for numerical fluxes.
 real(8), dimension(2) :: tdp
! Temporary storage for discontinuity postions.
 type(state) :: left_pressure, right_pressure
! A room for swap.
 type(state) :: left_room_for_swap, right_room_for_swap
end type comp_info

type flux_info
 type(state), dimension(4) :: l_flux, r_flux, flux
! The left and right fluxes in the critical cell.
 type(state), dimension(4,0:2) :: l_cff, l_cfb, r_cff, r_cfb
! The coefficients of the flux polynomials in the left and right
! sides.
 type(state), dimension(2) :: numerical_regularization
 real(8) :: numer_reg_coefficient
! Amount of numerical discrepancy on discontinuity curve and the numerical dissipation coefficient.
end type flux_info

public  clean_up_f_cell, clean_up_c_cell


contains


subroutine clean_up_f_cell(f_cell)

implicit none

type(flux_info), intent(out) :: f_cell
integer :: i, j

do i=1,4 
 f_cell%l_flux(i)=error_data; f_cell%r_flux(i)=error_data
 f_cell%flux(i)=error_data
 do j=0,2
  f_cell%l_cff(i,j)=error_data; f_cell%l_cfb(i,j)=error_data
  f_cell%r_cff(i,j)=error_data; f_cell%r_cfb(i,j)=error_data
 end do
end do

end subroutine clean_up_f_cell


subroutine clean_up_c_cell(c_cell)

implicit none

type(comp_info), intent(out) :: c_cell
integer :: i, j

do i=0,2
 c_cell%l_state(i)=error_data; c_cell%r_state(i)=error_data
 c_cell%or_state(i)=error_data
end do
do i=1,4
 do j=1,2
  c_cell%flux(i,j)=error_data
 end do
end do
c_cell%tdp=error_data

end subroutine clean_up_c_cell


end module computational_information