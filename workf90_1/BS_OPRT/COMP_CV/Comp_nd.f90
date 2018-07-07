module compute_node_ca
! This module computes cell-averages in node cells.

use node_cells
! nd_cls.f90'

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

use geo_info_cell, ng=>neighbor_edge
! 'gif_cell.f90'

use grid_map
! 'grid_map.f90'

implicit none

type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell
type(flux_info) :: f_cell
type(curve_plug), pointer :: tempp
type(cv_plug_info) :: plug
type(node_info) :: n_cell
type(state) :: uc, uc1
type(state), dimension(4) :: fl

private temp, g_cell, f_cell, tempp, plug, uc, uc1, &
        fl, n_cell
public  comput_ndca


contains


subroutine comput_ndca
! This subroutine computes the cell-average in a node cell

implicit none
integer :: i
type(state) :: fff

! call check_ring(1,'counter ')

do i=1,ndn
 if(ndd(i)%status.eq.'awake ') then
  call comput_ndfl2(ndd(i))
  uc=ndd(i)%n_cell%or_state
  fl=ndd(i)%n_cell%flux
  uc1=uc-r*(fl(2)-fl(4)+fl(3)-fl(1))
  ndd(i)%n_cell%or_state=uc1
  fff=ndd(i)%n_cell%or_state
 end if
end do

end subroutine comput_ndca


subroutine comput_ndfl2(ndd)
! This subroutine computes numerical fluxes across cell-edges 
! auxiliary critical cells.

implicit none
type(node_cell), intent(inout) :: ndd

integer :: head_mark, edge_c, edge_n
logical :: head_switch
character*3 :: empty_curve
integer :: i0, j0
type(adss_info) :: adss

n_cell=ndd%n_cell; fl=n_cell%flux; call clean_up_g_cell(g_cell)

tempp=>ndd%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.

do while(head_switch)
! Find the critical cell connected to the node cell.
 empty_curve='no '; plug=tempp%plug
 if(plug%end_type.eq.'begin ') then
  if(associated(acvv(plug%cv_nmb)%begin%next)) then
   temp=>acvv(plug%cv_nmb)%begin%next
   g_cell=temp%g_cell; edge_c=g_cell%edge(1)
  else
   empty_curve='yes'
  end if
 else
  if(associated(acvv(plug%cv_nmb)%eend%previous)) then
   temp=>acvv(plug%cv_nmb)%eend%previous
   g_cell=temp%g_cell; edge_c=g_cell%edge(2)
  else
   empty_curve='yes'
  end if
 end if

! Pick the numerical flux across the cell-edge.   
 if(empty_curve.eq.'no ') then
  g_cell=temp%g_cell; f_cell=temp%f_cell
  edge_n=neighbor_edge(edge_c)
  if(fl(edge_n).lt.0.9d0*error_data) then 
   fl(edge_n)=f_cell%flux(edge_c)   
  else
   fl(edge_n)=fl(edge_n)+f_cell%flux(edge_c) 
   if(tempp%next%address%idx.ne.head_mark) then
    if(tempp%plug%side_behind.eq.'left  ') then
     fl(edge_n)=fl(edge_n)-f_cell%l_flux(edge_c)
    else
     fl(edge_n)=fl(edge_n)-f_cell%r_flux(edge_c)
    end if
   else
	if(temp%l_stk.eq.adss_info(0,0).and.temp%r_stk.eq.adss_info(0,0)) then
	 call error_message; pause
    end if
	if(temp%l_stk.ne.adss_info(0,0)) then
	 fl(edge_n)=fl(edge_n)-f_cell%l_flux(edge_c)
    end if
    if(temp%r_stk.ne.adss_info(0,0)) then
	 fl(edge_n)=fl(edge_n)-f_cell%r_flux(edge_c)
    end if
   end if 		 		 	  
  end if
 else
  edge_n=plug%edge
  fl(edge_n)%value=0.d0
 end if 

 tempp=>tempp%next
 head_switch=(tempp%address%idx.ne.head_mark)

end do

! Compute the numerical flux across cell-edges adjacent to 
! critical cells that are not connected to the node cell. 
do edge_n=1, 4
 if(fl(edge_n).lt.0.9d0*error_data) then
! The flux across the cell-edge must not be computed yet.
  select case(edge_n)
   case(1); i0=n_cell%x_idx; j0=n_cell%y_idx-1; edge_c=3
   case(2); i0=n_cell%x_idx+1; j0=n_cell%y_idx; edge_c=4
   case(3); i0=n_cell%x_idx; j0=n_cell%y_idx+1; edge_c=1
   case(4); i0=n_cell%x_idx-1; j0=n_cell%y_idx; edge_c=2
  end select

  if(ggd_cell(i0,j0)%region.eq.'crit  ') then
   adss=ggd_cell(i0,j0)%ccpt(edge_n)%address
   call visit(adss,temp); f_cell=temp%f_cell
   fl(edge_n)=f_cell%flux(edge_c)
  else
   print*, 'There must be something wrong in node cell.'
   pause
  end if
 end if
end do    
 
n_cell%flux=fl; ndd%n_cell=n_cell

end subroutine comput_ndfl2

end module compute_node_ca