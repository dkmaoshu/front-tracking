module data_preparation
! This module prepares the smooth data used in the reset of
! discontinuity curves.

use solu_comput
! 'solu_com.f90'

use node_cells
! 'nd_cls.f90'

implicit none

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(phy_info) :: p_cell
!type(adss_info) :: adds_s

public  data_preparing
private data_prep, pick_stencil


contains 


subroutine data_preparing

implicit none
integer i

do i=1,cvn
 if(acvv(i)%status.eq.'awake') then
  call data_prep(acvv(i))
  call pick_stencil(acvv(i))
 end if
end do

end subroutine data_preparing


subroutine data_prep(acv)
! Prepare the local data of solution, 'u_exd', for the reset.

implicit none
type(auxi_discv), intent(inout) :: acv
! The single auxiliary discontinuity curve for which the data are to be prepared.

integer :: head_mark, end_mark, i, io, jo
logical :: head_switch
!type(crit_info_a) :: crit
type(state), dimension(-3:3) :: uu
type(state) :: fff
type(adss_info) :: ggg

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
 g_cell=temp%g_cell; p_cell=temp%p_cell
 fff=p_cell%l_state
 fff=P_cell%r_state
 ggg=temp%address

 do i=-1,1
  temp%ulx_exd(i)%value=error_data
  temp%uly_exd(i)%value=error_data
  temp%urx_exd(i)%value=error_data
  temp%ury_exd(i)%value=error_data
 end do
 temp%udiff_xl=error_data; temp%udiff_yl=error_data
 temp%udiff_xr=error_data; temp%udiff_yr=error_data
 temp%udiff_xl_f=error_data; temp%udiff_yl_f=error_data
 temp%udiff_xr_f=error_data; temp%udiff_yr_f=error_data
 temp%udiff_xl_b=error_data; temp%udiff_yl_b=error_data
 temp%udiff_xr_b=error_data; temp%udiff_yr_b=error_data


! The following computes the local smooth data associated to the critical cell.
 if(g_cell%g_type.eq.'xx '.or.g_cell%g_type.eq.'yy ') then
  call smth_dpr(io,jo,temp,uu,g_cell%g_type,'left  ','critical',.true.,3,'full    ')
  if(g_cell%g_type.eq.'xx ') temp%ulx_exd=uu(-1:1)
  if(g_cell%g_type.eq.'yy ') temp%uly_exd=uu(-1:1)
  call smth_dpr(io,jo,temp,uu,g_cell%g_type,'right ','critical',.true.,3,'full    ')
  if(g_cell%g_type.eq.'xx ') temp%urx_exd=uu(-1:1)
  if(g_cell%g_type.eq.'yy ') temp%ury_exd=uu(-1:1)
 else
  call smth_dpr(io,jo,temp,uu,'xx ','left  ','critical',.true.,3,'full    ')
  temp%ulx_exd=uu(-1:1)
  call smth_dpr(io,jo,temp,uu,'yy ','left  ','critical',.true.,3,'full    ')
  temp%uly_exd=uu(-1:1)
  call smth_dpr(io,jo,temp,uu,'xx ','right ','critical',.true.,3,'full    ')
  temp%urx_exd=uu(-1:1)
  call smth_dpr(io,jo,temp,uu,'yy ','right ','critical',.true.,3,'full    ')
  temp%ury_exd=uu(-1:1)
 end if

! The following computes the difference quotients in $x$- and $y$-directions.
 if(g_cell%g_type.eq.'xx '.or.g_cell%g_type.eq.'yy ') then
! Compute difference quotients only for critical cells of either 'xx'- or 'yy'-types.
  call smth_dpr(io,jo,temp,uu,'xx ','left  ','critical',.true.,2,'full    ')
  temp%udiff_xl=0.5d0*(uu(1)-uu(-1)) 
  temp%udiff_xl_f=uu(1)-uu(0); temp%udiff_xl_b=uu(0)-uu(-1) 
  call smth_dpr(io,jo,temp,uu,'yy ','left  ','critical',.true.,2,'full    ')
  temp%udiff_yl=0.5d0*(uu(1)-uu(-1)) 
  temp%udiff_yl_f=uu(1)-uu(0); temp%udiff_yl_b=uu(0)-uu(-1) 
  call smth_dpr(io,jo,temp,uu,'xx ','right ','critical',.true.,2,'full    ')
  temp%udiff_xr=0.5d0*(uu(1)-uu(-1))
  temp%udiff_xr_f=uu(1)-uu(0); temp%udiff_xr_b=uu(0)-uu(-1)    
  call smth_dpr(io,jo,temp,uu,'yy ','right ','critical',.true.,2,'full    ')
  temp%udiff_yr=0.5d0*(uu(1)-uu(-1)) 
  temp%udiff_yr_f=uu(1)-uu(0); temp%udiff_yr_b=uu(0)-uu(-1) 
 end if

 temp%tangled='no '

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if

end do

end subroutine data_prep


subroutine pick_stencil(acv)
! This subroutine picks the interpolation stencil for each critical cell.

implicit none
type(auxi_discv), intent(inout) :: acv
! The single auxiliary discontinuity curve to be updated.

integer :: head_mark, end_mark, i, k
logical :: head_switch
type(geo_info) :: gp_cell, gn_cell
real(8), dimension(0:3,2) :: temp_stencil
integer, dimension(0:3) :: temp_stencil_disp
real(8), dimension(3,2) :: stencil
integer, dimension(3) :: stencil_disp
real(8), dimension(2) :: pt1, pt2, pt, normal, normal_p, normal_n
real(8) :: difference_p, difference_n, length
character*8 :: chosen

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
 stencil=error_data; call clean_up_g_cell(gp_cell); call clean_up_g_cell(gn_cell)
 difference_p=error_data; difference_n=error_data
 temp_stencil=error_data; pt1=error_data; pt2=error_data; pt=error_data
 normal=error_data; normal_p=error_data; normal_n=error_data
 temp_stencil_disp=error_index; stencil_disp=error_index

 g_cell=temp%g_cell

! Pick stencils for critical cells.
  
 call pick_point(g_cell,1,pt1)
 temp_stencil(1,:)=pt1; temp_stencil_disp(1)=1
 call pick_point(g_cell,2,pt2)
 temp_stencil(2,:)=pt2; temp_stencil_disp(2)=2
! The two discontinuities in the critical cell under concern must be picked.

! The following picks the third point of the stencil. The third point exists only for 
! either 'xx'- or 'yy'-type critical cells. 
! First choose the critical cell.
 if(g_cell%g_type.eq.'xx '.or.g_cell%g_type.eq.'yy ') then
  if(associated(temp%p_nxt_dr)) gp_cell=temp%p_nxt_dr%g_cell
  if(associated(temp%n_nxt_dr)) gn_cell=temp%n_nxt_dr%g_cell

  call compute_normal(g_cell,normal)
  if(gp_cell%g_type.eq.'xx '.or.gp_cell%g_type.eq.'yy ') then
   call compute_normal(gp_cell,normal_p)
   difference_p=outer_product(normal_p,normal)
  end if
  if(gn_cell%g_type.eq.'xx '.or.gn_cell%g_type.eq.'yy ') then
   call compute_normal(gn_cell,normal_n)
   difference_n=outer_product(normal_n,normal)
  end if
  if(difference_p.ge.0.9d0*error_data) then
   if(difference_n.ge.0.9d0*error_data) then
    if(dabs(difference_p).lt.dabs(difference_n)) then
     chosen='previous'
    else
     chosen='next    '
    end if
   else
    chosen='previous'
   end if	  
  else
   if(difference_n.ge.0.9d0*error_data) then
    chosen='next    '
   else
    chosen='none    '
   end if
  end if   	 
! Then pick the point.
  select case(chosen)
   case('previous')
    call pick_point(gp_cell,1,pt)
    pt=pt+(/dfloat(gp_cell%x_idx-g_cell%x_idx),dfloat(gp_cell%y_idx-g_cell%y_idx)/)
    temp_stencil(0,:)=pt
   case('next    ')
    call pick_point(gn_cell,2,pt)
    pt=pt+(/dfloat(gn_cell%x_idx-g_cell%x_idx),dfloat(gn_cell%y_idx-g_cell%y_idx)/)
    temp_stencil(3,:)=pt
   case('none    ')
   case default; call error_message
  end select

 end if

! When the first two stencil points are two close, one of them must be deleted.
 length=length_of_segment(pt1,pt2)
 if(length.lt.1.0d-1) then
  select case(chosen)
   case('previous')
    temp_stencil(1,:)=temp_stencil(0,:); temp_stencil_disp(1)=error_index
	temp_stencil(0,:)=error_data; temp_stencil_disp(0)=error_index
   case('next    ')
    temp_stencil(2,:)=temp_stencil(3,:); temp_stencil_disp(2)=error_index
	temp_stencil(3,:)=error_data; temp_stencil_disp(3)=error_index
   case('none    ')
  end select
end if
  
! Finally, form the stencil.
 k=0
 do i=0,3
  if(temp_stencil(i,1).ge.0.9d0*error_data) then
   k=k+1
   stencil(k,:)=temp_stencil(i,:); stencil_disp(k)=temp_stencil_disp(i)
  end if
 end do
 if(k.gt.3) call error_message   
    
 temp%stencil=stencil; temp%stencil_disp=stencil_disp
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if

end do


contains


function outer_product(normal_1,normal_2) result (c)

implicit none
real(8), dimension(2), intent(in) :: normal_1, normal_2
real(8) :: c

c=normal_1(1)*normal_2(2)-normal_1(2)*normal_2(1)

end function outer_product


subroutine compute_normal(g_cell,normal)

implicit none
type(geo_info), intent(in) :: g_cell
real(8), dimension(2), intent(out) :: normal

real(8), dimension(2) :: pt1, pt2

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)
normal=normal_of_line(pt1,pt2)

end subroutine compute_normal


end subroutine pick_stencil


end module data_preparation