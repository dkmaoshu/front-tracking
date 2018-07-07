module compute_cell_average
! This module is for the computation of ordinary cell-averages
! in critical cells.

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

use geo_info_cell  !, n_g=>neighbor_edge
! gif_cell.f90'

implicit none

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell, ng_cell
type(phy_info) :: p_cell, np_cell
type(flux_info) :: f_cell, nf_cell
type(comp_info) :: c_cell, nc_cell
type(state) :: uc, uc1
type(state), dimension(4) :: fl

private temp, g_cell, ng_cell, p_cell, f_cell, comp_flux, &
        flux_int, uc, uc1, fl, comput_cvca, np_cell 
public  compute_cvcas


contains


subroutine compute_cvcas
! This subroutine computes cell-avarages on auxiliary discontinuity
! curves. The parameter 'level' stands for the stage in the
! Predictor-corrector procedure.

implicit none

integer :: i

! call check_list_ac(1,'down  ')

do i=1, cvn
 if(acvv(i)%status.eq.'awake') then
  call comput_cvca(acvv(i))
 end if
end do

end subroutine compute_cvcas


subroutine comput_cvca(acv)
! Computation on one auxiliary curve.

implicit none

type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

! type(state) :: fff, ggg, hhh

nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
 boundary_end=>acv%eend_boundary%previous
 if(acv%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acv%eend%previous%next%address%idx
 end if
end if

do while(associated(temp).and.head_switch) 
 address=temp%address
 p_cell=temp%p_cell; g_cell=temp%g_cell
!  fff=p_cell%l_state; ggg=p_cell%r_state
 f_cell=temp%f_cell; c_cell=temp%c_cell
 uc=p_cell%or_state
 
 call clean_up_g_cell(ng_cell); call clean_up_p_cell(np_cell)

!  call find_state(p_cell%l_state,fff)

 call comp_flux
 uc1=uc-r*(fl(2)-fl(4)+fl(3)-fl(1))
!  hhh=uc1
  
!  fff=fl(2)-fl(4)+fl(3)-fl(1)

 temp%p_cell%or_state=uc1
 temp%f_cell%flux=fl

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine comput_cvca


subroutine comp_flux
! Compute the numerical fluxes across the four cell-edges of an
! auxiliary critical cell.
implicit none

integer :: edge, neighbor_case, i, j, point, edge_comput
character*6 :: location
type(state) :: xx1, xx2
real(8) :: position
character*8 :: partner

 type(state) :: xxxxx, yyyyy

do i=1,4; fl(i)=error_data; end do

! The following computes the numerical fluxes across the cell-
! edges not crossed by the discontinuity curve.
do i=1,4
 j=cycle_4(i+1)
 if(g_cell%point(i).eq.g_cell%point(j)) then
  if(g_cell%point(i).eq.'left') then
   fl(i)=f_cell%l_flux(i)
  else
   if(g_cell%point(i).eq.'right') then
    fl(i)=f_cell%r_flux(i)
   else
    print*, 'There must be something wrong in comp_flux.'
    stop
   end if
  end if
 end if
end do

! The following computes the numerical fluxes across the cell
! edges crossed by the discontinuity curve.
partner=temp%partner
do i=1,2
 call clean_up_g_cell(ng_cell); call clean_up_f_cell(nf_cell)
 call clean_up_p_cell(np_cell); call clean_up_c_cell(nc_cell)
 edge=g_cell%edge(i); neighbor_case=1
 location='middle'; 
 position=error_data
 if(i.eq.1) then
  if(associated(temp%p_nxt_dr)) then
   ng_cell=temp%p_nxt_dr%g_cell; nf_cell=temp%p_nxt_dr%f_cell
   np_cell=temp%p_nxt_dr%p_cell; nc_cell=temp%p_nxt_dr%c_cell
   edge_comput=neighbor_edge(edge)
  else
!   location='end   '
   if(partner.eq.'previous') then
    ng_cell=temp%previous%g_cell; nf_cell=temp%previous%f_cell
	np_cell=temp%previous%p_cell; nc_cell=temp%previous%c_cell
	edge_comput=edge
   else
    location='end   '
   end if
  end if
 else
  if(associated(temp%n_nxt_dr)) then
   ng_cell=temp%n_nxt_dr%g_cell; nf_cell=temp%n_nxt_dr%f_cell
   np_cell=temp%n_nxt_dr%p_cell; nc_cell=temp%n_nxt_dr%c_cell
   edge_comput=neighbor_edge(edge)
  else
!   location='end   '
   if(partner.eq.'next    ') then
    ng_cell=temp%next%g_cell; nf_cell=temp%next%f_cell
	np_cell=temp%next%p_cell; nc_cell=temp%next%c_cell
	edge_comput=edge
   else
    location='end   '
   end if
  end if
 end if
! Pick the next-door auxiliary critical cell.

 if(location.eq.'middle') then
  select case(g_cell%g_type) 
   case('xx')
    if(ng_cell%x_idx.gt.g_cell%x_idx) neighbor_case=2
    if(ng_cell%x_idx.lt.g_cell%x_idx) neighbor_case=3
   case('yy')
    if(ng_cell%y_idx.gt.g_cell%y_idx) neighbor_case=2
    if(ng_cell%y_idx.lt.g_cell%y_idx) neighbor_case=3
  end select
 end if
! Determine the neighboring case, the way the two auxiliary 
! critical cells are next-door neighbored.

 position=0.5d0*(g_cell%dis_posi(i)+c_cell%tdp(i))
 select case(neighbor_case)
  case(2); position=position-1.0d0
  case(3); position=position+1.0d0
 end select
! Compute 'position' according to the wind direction.

 select case(edge)
  case(1,2); point=edge
  case(3,4); point=cycle_4(edge+1)
 end select

 select case(neighbor_case)

  case(1) 

  if(g_cell%point(point).eq.'left') then
   xx1=flux_int(f_cell%l_cfb,position,edge)- &
   flux_int(f_cell%l_cfb,-0.5d0,edge)
   xx2=flux_int(f_cell%r_cff,0.5d0,edge)- &
   flux_int(f_cell%r_cff,position,edge)
     
    xxxxx=flux_int(f_cell%r_cff,0.5d0,edge); yyyyy=flux_int(f_cell%r_cff,position,edge)
    
  else
   xx1=flux_int(f_cell%r_cfb,position,edge)- &
   flux_int(f_cell%r_cfb,-0.5d0,edge)
   xx2=flux_int(f_cell%l_cff,0.5d0,edge)- &
   flux_int(f_cell%l_cff,position,edge)
    
    xxxxx=flux_int(f_cell%l_cff,0.5d0,edge); yyyyy=flux_int(f_cell%l_cff,position,edge)
    
  end if
   
 case(2) 

  if(g_cell%point(point).eq.'left') then
   xx1=flux_int(nf_cell%l_cfb,position,edge_comput)- &
   flux_int(nf_cell%l_cfb,-0.5d0,edge_comput)
   xx1=xx1+f_cell%l_flux(edge)
   xx2=flux_int(nf_cell%r_cff,-0.5d0,edge_comput)- &
   flux_int(nf_cell%r_cff,position,edge_comput)
  else
   xx1=flux_int(nf_cell%r_cfb,position,edge_comput)- &
   flux_int(nf_cell%r_cfb,-0.5d0,edge_comput)
   xx1=xx1+f_cell%r_flux(edge)
   xx2=flux_int(nf_cell%l_cff,-0.5d0,edge_comput)- &
   flux_int(nf_cell%l_cff,position,edge_comput)
  end if
  
  case(3)

  if(g_cell%point(point).eq.'left') then
   xx1=flux_int(nf_cell%l_cfb,position,edge_comput)- &
   flux_int(nf_cell%l_cfb,0.5d0,edge_comput)
   xx2=flux_int(nf_cell%r_cff,0.5d0,edge_comput)- &
   flux_int(nf_cell%r_cff,position,edge_comput)
   xx2=xx2+f_cell%r_flux(edge)
  else
   xx1=flux_int(nf_cell%r_cfb,position,edge_comput)- &
   flux_int(nf_cell%r_cfb,0.5d0,edge_comput)
   xx2=flux_int(nf_cell%l_cff,0.5d0,edge_comput)- &
   flux_int(nf_cell%l_cff,position,edge_comput)
   xx2=xx2+f_cell%l_flux(edge)
  end if
  
 end select

 fl(edge)=xx1+xx2
 ! Compute the flux for each case.
end do

end subroutine comp_flux


function flux_int(flux_cf,x,edge_index)
! The function to compute the flux coefficients used in the 
! computation of the numerical fluxes across the cell-edges 
! crossed by the discontinuity curve.

implicit none
type(state), intent(in) :: flux_cf(4,0:2)
type(state) :: flux_int
real(8), intent(in) :: x
integer, intent(in) :: edge_index

flux_int=flux_cf(edge_index,0)+x*flux_cf(edge_index,1)+ &
x*x*flux_cf(edge_index,2)

end function flux_int


end module compute_cell_average