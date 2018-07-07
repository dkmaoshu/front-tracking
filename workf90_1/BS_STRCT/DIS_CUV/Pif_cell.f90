module phy_info_cell
! This module describes the data structure of the physical
! information contained in critical cells.

use physical_state
implicit none

type phy_info
 
 type(state) :: l_state, r_state, or_state
 integer :: wv_nb
! Physical information consists of the left state 'l_state', 
! right state 'r_state', ordinary state 'or_state', and the 
! wave number 'wv_nb'.

end type phy_info

interface assignment(=)
 module procedure assignn
end interface

private assignn
public  clean_up_p_cell


contains


subroutine clean_up_p_cell(p_cell)

implicit none
type(phy_info) :: p_cell

p_cell%l_state%value=error_data
p_cell%r_state%value=error_data
p_cell%or_state%value=error_data
p_cell%wv_nb=-1000

end subroutine clean_up_p_cell


subroutine assignn(p_cell2,p_cell1)

implicit none
type(phy_info), intent(in) :: p_cell1
type(phy_info), intent(out) :: p_cell2

p_cell2%l_state=p_cell1%l_state
p_cell2%r_state=p_cell1%r_state
p_cell2%or_state=p_cell1%or_state
p_cell2%wv_nb=p_cell1%wv_nb

end subroutine assignn


end module phy_info_cell