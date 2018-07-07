module upd_cv_cases

use solu_comput
! 'solu_com.f90'

use discontinuity_curves
! 'discontinuity_curves.f90'

use tools_for_reset
! 'tools_rt.f90'

use numerical_integral_in_2D
! 'num_int2.f90'

implicit none

type(state), dimension(-1:1) :: ulx_exd_p, uly_exd_p, &
                                urx_exd_p, ury_exd_p

public  move, split, set_move_1, ulx_exd_p, uly_exd_p, &
      	urx_exd_p, ury_exd_p
private	geo_move, geo_split, integral_2, integral_1, reduced_function_in_1D, &
        function_in_2D, upper_limit, lower_limit, u_constant, slop_x, slop_y, &
		a_upper, b_upper, c_upper, a_lower, b_lower, c_lower, geometrical_type, &
		upper, lower, xy


contains


subroutine move(temp,g_cell,p_cell,uu)
! This subroutine moves auxiliary critical cell both geometrically 
! and physically. The updated information is stored in 'g_cell' 
! and 'p_cell', and the date for smooth solution in original grid 
! cell is stored in 'uu'. The update of neighbor relation and the
! grid-map is not concerned.

implicit none

type(auxi_crit_cell), pointer :: temp
! Pointing to the critical cell under concern.
type(geo_info), intent(out) :: g_cell
! Geometrical information of the moved critical cell.
type(phy_info), intent(out) :: p_cell
! Physical information of the moved critical cell.
type(state), intent(out) :: uu
! The state left behind by the move in the original grid cell.

type(auxi_crit_cell), pointer :: tempp
type(geo_info) :: g_tp
type(phy_info) :: p_tp
integer :: single, io, jo, in, jn, ij
real(8), dimension(2) :: pt1, pt2, normal, nmp, nmn
!real(8) :: velo
type(state) :: uu_l, uu_r
type(state) :: uulp, uurp, uuln, uurn, uul, uur
!character*12 :: wave_type

! type(geo_info) :: g_pp

! Initialization.
g_cell=temp%g_cell; p_cell=temp%p_cell; g_tp=g_cell
io=g_cell%x_idx; jo=g_cell%y_idx
call geo_move(g_tp)
in=g_tp%x_idx; jn=g_tp%y_idx

select case(g_cell%g_type)

 case('xx ','yy ')
! Compute the left and right states for the moved critical cell
! of 'xx' or 'yy' type.
  call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2)
  normal=normal_of_line(pt1,pt2)
  ij=in-io+jn-jo
  if(io.ne.in) then
   uu_l=temp%ulx_exd(ij); uu_r=temp%urx_exd(ij)
  else 
   uu_l=temp%uly_exd(ij); uu_r=temp%ury_exd(ij)
  end if
!  call riemann(uu_l,uu_r,normal,p_cell%wv_nb,uul,uur, &
!               velo,wave_type)
  uul=uu_l; uur=uu_r
 
 case('xy ')
! Compute the left and right states for the moved critical cell
! of 'xy' type.  
  nmp=error_data; nmn=error_data
  uulp%value=error_data; uurp%value=error_data
  uuln%value=error_data; uurn%value=error_data

  call find_single(g_cell,single)
  select case(single)
   case(1); in=io-1; jn=jo-1
   case(2); in=io+1; jn=jo-1
   case(3); in=io+1; jn=jo+1
   case(4); in=io-1; jn=jo+1
  end select

! For the moved critical cell of 'xy' type, extrapolation in 
! two directions should be carried out to compute the left
! and right states.      
  if(associated(temp%next).and.temp%next%g_cell%g_type &
     .ne.'xy ') then
   tempp=>temp%next
   call pick_point(g_cell,1,pt1)
   call pick_point(tempp%g_cell,2,pt2)
   nmn=normal_of_line(pt1,pt2)
   if(tempp%g_cell%g_type.eq.'xx ') then
    ij=in-io
   else
    ij=jn-jo
   end if
   if(io.eq.tempp%g_cell%x_idx) then
    uu_l=tempp%ulx_exd(ij); uu_r=tempp%urx_exd(ij)
   else
    uu_l=tempp%uly_exd(ij); uu_r=tempp%ury_exd(ij)
   end if    
!   call riemann(uu_l,uu_r,nmn,p_cell%wv_nb,uuln,uurn, &
!                velo,wave_type)
   uuln=uu_l; uurn=uu_r
  end if 

  if(associated(temp%previous)) then
   tempp=>temp%previous
   call pick_point(g_cell,2,pt2)
   call pick_point(tempp%g_cell,1,pt1)
   nmp=normal_of_line(pt1,pt2)
   if(tempp%g_cell%g_type.eq.'xx ') then
    ij=in-io
   else
    ij=jn-jo
   end if
   if(in.eq.tempp%g_cell%x_idx) then
    uu_l=ulx_exd_p(ij); uu_r=urx_exd_p(ij)
   else
    uu_l=uly_exd_p(ij); uu_r=ury_exd_p(ij)
   end if
   if(uu_l.gt.0.9d0*error_data.and.uu_r.gt.0.9d0*error_data) then
!    call riemann(uu_l,uu_r,nmp,p_cell%wv_nb,uulp,uurp, &
!                 velo,wave_type)
   uulp=uu_l; uurp=uu_r
   end if
  end if

! The final left and right states should be taken as the
! averages of the results obtained in the two directions.
  if(uulp.gt.0.9d0*error_data) then
   uul=uulp; uur=uurp
   if(uuln.gt.0.9d0*error_data) then
    uul=0.5d0*(uul+uuln); uur=0.5d0*(uur+uurn)
   end if											    
  else
   if(uuln.gt.0.5d0*error_data) then
    uul=uuln; uur=uurn
   else
    uul=p_cell%l_state; uur=p_cell%r_state
   end if
  end if
     
 end select

! Compute the left, right and ordinary states.
 p_tp%l_state=uul; p_tp%r_state=uur; p_tp%wv_nb=p_cell%wv_nb
 if(move_dir(1).eq.'l_to_r') then
  p_tp%or_state=p_cell%or_state-p_cell%l_state+uur
  uu=p_cell%l_state
 else
  p_tp%or_state=p_cell%or_state-p_cell%r_state+uul
  uu=p_cell%r_state
 end if							  
 
 g_cell=g_tp; p_cell=p_tp

end subroutine move


subroutine split(temp,g_cell,g_new,p_cell,p_new,num)
! This subroutine splits auxiliary critical cell both geometrically
! and physically. The old and new geometrical and physical 
! information is stored in 'g_cell', g_new', 'p_cell' and 'p_new'.
! 'num' is the number of the discontinuity position causing the
! split. The update of neighbor relation and grid-map is not
! concerned.

implicit none

type(auxi_crit_cell), pointer :: temp
! Pointing to the critical cell under concern.
type(geo_info), intent(out) :: g_cell, g_new 
! Geometric information of the split and new critical cell 
! produced in the split.
type(phy_info), intent(out) :: p_cell, p_new
! Physical infomation of the split and new critical cell
! produced in the split.
integer, intent(in) :: num
! The number of the discontinuity position causing the split.

integer :: io, jo, in, jn, ij, i, edge_num
integer, dimension(2) :: edge
real(8), dimension(2) :: normal, coeff
real(8), dimension(3,2) :: pt1, pt2
real(8) :: c_between !, velo
type(state) :: uu_l, uu_r !, uul, uur
type(state), dimension(2) :: difference, uu_old, uu_new
!character*12 :: wave_type
real(8), dimension(3,2) :: stencil
real(8), dimension(3) :: cv_coeffs
character*3 :: g_type_mem

! Initialization.
call clean_up_g_cell(g_new); call clean_up_p_cell(p_new)
g_cell=temp%g_cell; p_cell=temp%p_cell; stencil=temp%stencil
!geometrical_type=g_cell%g_type
edge=g_cell%edge; edge_num=g_cell%edge(num); g_type_mem=g_cell%g_type

! First determine the geometrical_type.
select case(edge_num)
 case(1,3); geometrical_type='xx '
 case(2,4); geometrical_type='yy '
 case default; call error_message
end select

! Compute the normal of the discontinuity segment.
call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2)
normal=normal_of_line(pt1,pt2)

! Computed the area of portion of the critical cell that splits 
! off.
call geo_split(stencil,g_cell,g_new,num,cv_coeffs)

! The following determine the 2D integrant.
u_constant=p_cell%r_state-p_cell%l_state
if(g_cell%x_idx.eq.g_new%x_idx) then
 if(g_cell%y_idx.eq.g_new%y_idx) call error_message
 slop_x=temp%udiff_xr-temp%udiff_xl
 if(g_cell%y_idx.gt.g_new%y_idx) then
  slop_y=temp%udiff_yr_b-temp%udiff_yl_b
 else
  slop_y=temp%udiff_yr_f-temp%udiff_yl_f
 end if    
else
 if(g_cell%y_idx.ne.g_new%y_idx) call error_message 
 slop_y=temp%udiff_yr-temp%udiff_yl
 if(g_cell%x_idx.gt.g_new%x_idx) then
  slop_x=temp%udiff_xr_b-temp%udiff_xl_b
 else
  slop_x=temp%udiff_xr_f-temp%udiff_xl_f
 end if
end if

do i=1,2
 a_upper=error_data; b_upper=error_data; c_upper=error_data
 a_lower=error_data; b_lower=error_data; c_lower=error_data

! The following determines the upper and lower limits of the integral.
 if(i.eq.num) then 
  select case(geometrical_type)
   case('xx ')
    c_between=0.5d0*dfloat(g_new%x_idx-g_cell%x_idx)
    if(g_new%x_idx.gt.g_cell%x_idx) then
!     c_between=0.5d0*dfloat(g_new%x_idx-g_cell%x_idx)
     call set_limits('curve_top   ')
    else
!     c_between=0.5d0*dfloat(g_cell%x_idx-g_new%x_idx)
     call set_limits('curve_below ')
    end if
   case('yy ')
    c_between=0.5d0*dfloat(g_new%y_idx-g_cell%y_idx)
    if(g_new%y_idx.gt.g_cell%y_idx) then
!     c_between=0.5d0*dfloat(g_new%y_idx-g_cell%y_idx)
     call set_limits('curve_top   ')
    else
!     c_between=0.5d0*dfloat(g_cell%y_idx-g_new%y_idx)
     call set_limits('curve_below ')
    end if
   case default; call error_message
  end select
 else
  select case(geometrical_type)
   case('xx ')
    c_between=0.5d0*dfloat(g_new%x_idx-g_cell%x_idx)
    if(g_new%x_idx.gt.g_cell%x_idx) then
!     c_between=0.5d0*dfloat(g_new%x_idx-g_cell%x_idx)
     call set_limits('curve_below ')
    else
!     c_between=0.5d0*dfloat(g_cell%x_idx-g_new%x_idx)
     call set_limits('curve_top   ')
    end if
   case('yy ')
    c_between=0.5d0*dfloat(g_new%y_idx-g_cell%y_idx)
    if(g_new%y_idx.gt.g_cell%y_idx) then
!     c_between=0.5d0*dfloat(g_new%y_idx-g_cell%y_idx)
     call set_limits('curve_below ')
    else
!     c_between=0.5d0*dfloat(g_cell%y_idx-g_new%y_idx)
     call set_limits('curve_top   ')
    end if
   case default; call error_message
  end select
 end if

! The following determines the upper and lower.
if(edge(i).eq.edge_num) then
 select case(edge_num)
  case(2,3) 
   lower=g_cell%dis_posi(num); upper=0.5d0
  case(1,4)
   lower=-0.5d0; upper=g_cell%dis_posi(num)
 end select
else
 select case(edge_num)
  case(2,3)
   lower=-0.5d0; upper=g_cell%dis_posi(num)
  case(1,4) 
   lower=g_cell%dis_posi(num); upper=0.5d0
 end select
end if

! The following computes the difference.
 call integral_2(difference(i))
end do

! Find the left and right states for the new critical cell.
io=g_cell%x_idx; jo=g_cell%y_idx
in=g_new%x_idx; jn=g_new%y_idx
ij=in-io+jn-jo
if(io.ne.in) then
 uu_l=temp%ulx_exd(ij); uu_r=temp%urx_exd(ij)
else
 uu_l=temp%uly_exd(ij); uu_r=temp%ury_exd(ij)
end if 
!call riemann(uu_l,uu_r,normal,p_cell%wv_nb,uul,uur, &
!             velo,wave_type)
!p_new%l_state=uul; p_new%r_state=uur; p_new%wv_nb=p_cell%wv_nb
p_new%l_state=uu_l; p_new%r_state=uu_r; p_new%wv_nb=p_cell%wv_nb

! Compute the 'or_state' for the new critical cell and update 
! the 'or_state' for the old one.
select case(move_dir(num))
 case('l_to_r')
  uu_old(num)=p_cell%or_state+difference(num)
  uu_new(num)=p_new%r_state-difference(num)
  uu_old(3-num)=p_cell%l_state+difference(3-num)
  uu_new(3-num)=p_cell%or_state+p_new%r_state-p_cell%l_state-difference(3-num)
 case('r_to_l')
  uu_old(num)=p_cell%or_state-difference(num)
  uu_new(num)=p_new%l_state+difference(num)
  uu_old(3-num)=p_cell%r_state-difference(3-num)
  uu_new(3-num)=p_cell%or_state+p_new%l_state-p_cell%r_state+difference(3-num)
 case default; call error_message
end select


select case(edge_num)
 case(2,3)
  coeff(num)=0.5d0-g_cell%dis_posi(num)
  coeff(3-num)=g_cell%dis_posi(num)+0.5d0
 case(1,4)
  coeff(num)=g_cell%dis_posi(num)+0.5d0
  coeff(3-num)=0.5d0-g_cell%dis_posi(num)
end select
  
!do i=1,2
! select case(edge(i))
!  case(2,3)
!   coeff(i)=0.5d0-g_cell%dis_posi(num)
!   coeff(3-i)=g_cell%dis_posi(num)+0.5d0
!  case(1,4)
!   coeff(i)=g_cell%dis_posi(num)+0.5d0
!   coeff(3-i)=0.5d0-g_cell%dis_posi(num)
! end select
!end do

p_cell%or_state=coeff(2)*uu_old(1)+coeff(1)*uu_old(2)
p_new%or_state=coeff(2)*uu_new(1)+coeff(1)*uu_new(2)

! Special treatment for 'xy'-type critical cell.
if(g_type_mem.eq.'xy ') then
 p_cell%or_state=uu_old(num)
 p_new%or_state=uu_new(num)
end if

contains


subroutine set_limits(situation)

implicit none
character*12, intent(in) :: situation

select case(situation)
 case('curve_top   ')
  a_upper=cv_coeffs(3); b_upper=cv_coeffs(2); c_upper=cv_coeffs(1)
  a_lower=0.0d0; b_lower=0.0d0; c_lower=c_between
 case('curve_below ')  
  a_lower=cv_coeffs(3); b_lower=cv_coeffs(2); c_lower=cv_coeffs(1)
  a_upper=0.0d0; b_upper=0.0d0; c_upper=c_between
 case default; call error_message
end select

end subroutine set_limits


end subroutine split


end module upd_cv_cases