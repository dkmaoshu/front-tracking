module auxi_cc_prod
! This module describes the production of an auxiliary critical
! cell from two neighbouring critical cells.

use	geo_info_cell
! 'gif_cell.f90'

use phy_info_cell
! 'pif_cell.f90'

public with_geo, with_phy
private select_st


contains


subroutine with_geo(g_cell,gn_cell,ag_cell,neighbour)
! To produce the geometry cell of the auxiliary critical cell.

implicit none

type(geo_info) :: g_cell, gn_cell, ag_cell
! 'g_cell' is the geometry cell of the current critical cell,
! 'gn_cell' is the geometry cell of the neighbour critical cell, 
! and 'ag_cell' is the geometry cell of the auxiliary critical 
! cell to be produced.
character*8 :: neighbour
! To indicate whether the neighbour critical cell is the previous
! or the next one; the possible values for this variable is
! 'previous' and 'next'.
intent(in) :: g_cell, gn_cell
intent(out) :: ag_cell

integer :: first, second
! Variables 'first' and 'second' are the indexes of the two
! different edges or discontinuity positions in the two
! neighbouring critical cells. The 'first' one is in the 
! current critical cell and the 'second is in the neighbour
! critical cell.

call clean_up_g_cell(ag_cell)
! When no auxiliary critical cell is produced by the two neighbour
! critical cells, the geometry cell of the auxiliary critical
! cell takes the value set in the subroutine 'clean_up_g_cell'
if(gn_cell%g_type=='xy') then 

 if(neighbour.eq.'previous') then
  first=2; second=1
 else 
  first=1; second=2
 end if

 if(g_cell%edge(first).ne.gn_cell%edge(second)) then	
  
  ag_cell=g_cell

  ag_cell%edge(second)=gn_cell%edge(second)
  if(g_cell%x_idx==gn_cell%x_idx) then
   ag_cell%g_type='yy'       
   if(g_cell%y_idx>gn_cell%y_idx) then
    ag_cell%point(1)=side_shift(g_cell%point(4))
    ag_cell%point(2)=side_shift(g_cell%point(4))
    ag_cell%dis_posi(second)=gn_cell%dis_posi(second)-1.0d0
   else								  		          
    ag_cell%point(3)=side_shift(g_cell%point(1)) 
    ag_cell%point(4)=side_shift(g_cell%point(1))
    ag_cell%dis_posi(second)=gn_cell%dis_posi(second)+1.0d0
   end if
  else
   ag_cell%g_type='xx'       
   if(g_cell%x_idx>gn_cell%x_idx) then
    ag_cell%point(1)=side_shift(g_cell%point(3))
    ag_cell%point(4)=side_shift(g_cell%point(3))
    ag_cell%dis_posi(second)=gn_cell%dis_posi(second)-1.0d0
   else								  		          
    ag_cell%point(2)=side_shift(g_cell%point(1)) 
    ag_cell%point(3)=side_shift(g_cell%point(1))
    ag_cell%dis_posi(second)=gn_cell%dis_posi(second)+1.0d0
   end if
  end if      	
 end if
        	
end if
  
end subroutine with_geo


subroutine with_phy(g_cell,gn_cell,p_cell,pn_cell,ap_cell)
! To produce the physics cell of the auxiliary critical cell.

implicit none

type(geo_info) :: g_cell, gn_cell
! g_cell is the geometry cell of the current critical cell and
! gn_cell is the geometry cell of the neighbour critical cell.
type(phy_info) :: p_cell, pn_cell, ap_cell
! p_cell is the physics cell of the current critical cell,
! gp_cell is the physics cell of the neighbour critical cell, 
! and ap_cell is the physics cell of the auxiliary critical 
! cell to be produced.
intent(in) :: p_cell, pn_cell, g_cell, gn_cell
intent(out) :: ap_cell

ap_cell=p_cell
ap_cell%or_state=ap_cell%or_state+pn_cell%or_state
if(g_cell%x_idx.gt.gn_cell%x_idx.or.g_cell%y_idx.gt. &
gn_cell%y_idx) then
 ap_cell%or_state= &
 ap_cell%or_state-select_st(side_shift(g_cell%point(3)),pn_cell)
else
 ap_cell%or_state= &
 ap_cell%or_state-select_st(side_shift(g_cell%point(1)),pn_cell)
end if     

end subroutine with_phy


function select_st(side,p_cell)
implicit none
type(phy_info) :: p_cell
character*6 :: side
type(state) :: select_st

if(side.ne.'left'.and.side.ne.'right') then
 print*, 'There must be something wrong with the point-side'
 pause
end if
if(side.eq.'left') then
 select_st=p_cell%l_state
else
 select_st=p_cell%r_state
end if
  
end function select_st


end module auxi_cc_prod