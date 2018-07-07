module compute_discontinuity_positions
! This module is for the computation of discontinuity positions
! form middle discontinuity positions.

use pick_middle_positions
! 'pick_mdp.f90'

implicit none

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell

public  compute_positions
private temp, tempa, g_cell, pick_positions, compute_posi 


contains


subroutine compute_positions
! This subroutine executes the computation on all discontinuity
! curves.

implicit none

integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake') then
  call compute_posi(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

do i=1, cvn
 if(acvv(i)%status.eq.'awake') then
  call compute_xy_mid_posi(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

end subroutine compute_positions


subroutine compute_posi(acv)
! This subroutine preformances the computation on a single
! discontinuity curve.

implicit none
type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
character*8 :: partner
integer :: head_mark, end_mark
logical :: head_switch

! type(auxi_crit_cell), pointer :: temp_g
! type(geo_info) :: g_show

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
 g_cell=temp%g_cell; partner=temp%partner
 
 tempa=>temp
 call compute_position

 temp%g_cell=g_cell

! call visit(adss_info(1,13),temp_g)
! g_show=temp_g%g_cell

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do


contains


subroutine compute_position
! This subroutine preforms the computation in a critical cell.
! As a matter of fact, it computes only the second discontinuity
! position; the first position has already been computed in the
! previous next-door critical cell as the second position there.

implicit none

type(auxi_crit_cell), pointer :: tempap
type(geo_info) :: gp_cell
integer :: i

call clean_up_g_cell(gp_cell)

select case(partner)

 case('previous')
! Transfer the discontinuity positions from the previous partner.
  tempap=>tempa%previous; gp_cell=tempap%g_cell
  g_cell%dis_posi=gp_cell%dis_posi; g_cell%normal=gp_cell%normal
  do i=1,2
   select case(g_cell%edge(i))
    case(1,3)
     g_cell%dis_posi(i)=g_cell%dis_posi(i)+dfloat(gp_cell%x_idx-g_cell%x_idx)
    case(2,4)
     g_cell%dis_posi(i)=g_cell%dis_posi(i)+dfloat(gp_cell%y_idx-g_cell%y_idx)
   end select
  end do
  tempa%g_cell=g_cell
 
 case('next    ','single  ')
  if(.not.associated(temp%p_nxt_dr)) call compute(1)

  call compute(2)
  if(associated(tempa%n_nxt_dr)) then
   tempa%g_cell=g_cell
   tempap=>tempa%n_nxt_dr; gp_cell=tempap%g_cell
   gp_cell%dis_posi(1)=g_cell%dis_posi(2)
   gp_cell%normal(1,:)=g_cell%normal(2,:)
   select case(gp_cell%edge(1))
    case(1,3)
     gp_cell%dis_posi(1)=gp_cell%dis_posi(1)+dfloat(g_cell%x_idx-gp_cell%x_idx)
    case(2,4)
     gp_cell%dis_posi(1)=gp_cell%dis_posi(1)+dfloat(g_cell%y_idx-gp_cell%y_idx)
   end select
   tempap%g_cell=gp_cell
  end if
end select


end subroutine compute_position


subroutine compute(num)

implicit none
integer, intent(in) :: num
! The number of the edge of discontinuity.

real(8), dimension(2,2) :: positions
integer :: point_number
real(8), dimension(2) :: stencil, points, coefficients, &
                         pt1, pt2, normal
character*8 :: fashion
real(8) :: dis_posi, paramet
type(geo_info) :: gn
character*12 :: up_down_stream

! Pick the extrapolation/interpolation stencils and so on.
select case(num)
 case(1)
  call pick_positions('up    ',positions,point_number,fashion,up_down_stream)
 case(2)
  call pick_positions('down  ',positions,point_number,fashion,up_down_stream)
end select

! if(num.eq.1) print*, 'I have passed this point!'

select case(point_number)

 case(1)
! Stencil of only one point, must be a critical cell of 'xy'-type.
  if(g_cell%g_type.ne.'xy ') call error_message
  select case(fashion)
   case('single  ')
    call normal_of_single(g_cell,normal)
    g_cell%normal(num,:)=normal
    g_cell%dis_posi(num)=g_cell%mdis_posi(num)
   case('corner  ')
	call normal_of_corner(num,normal)
	g_cell%normal(num,:)=normal
	g_cell%dis_posi(num)=g_cell%mdis_posi(num)
   case default; call error_message
  end select

 case(2)
! Stencil of two points.  
! First compute the normal.
  select case(up_down_stream)
   case('down_stream '); pt1=positions(1,:); pt2=positions(2,:)
   case('up_stream   '); pt1=positions(2,:); pt2=positions(1,:)
   case default; call error_message 
  end select

  normal=normal_of_line(pt1,pt2)
  g_cell%normal(num,:)=normal 

! Then compute the discontinuity postions.
  select case(g_cell%edge(num))
   case(1,3)
    points=positions(:,1)
    stencil=positions(:,2)
   case(2,4)
    stencil=positions(:,1)
    points=positions(:,2)
  end select

! Form the polynomials corresponding to the x- or y-directions.
  call polynomial(points,stencil,coefficients)

! Compute the discontinuity position.
  select case(g_cell%edge(num))
   case(1,4); paramet=-0.5d0
   case(2,3); paramet=0.5d0  
  end select
  dis_posi=polynomial_function(coefficients,paramet)
  g_cell%dis_posi(num)=dis_posi

end select

if(g_cell%g_type.eq.'xy ') g_cell%dis_posi(num)=g_cell%mdis_posi(num)
if(num.eq.2) then
 if(associated(temp%n_nxt_dr)) then
  gn=temp%n_nxt_dr%g_cell
  if(gn%g_type.eq.'xy ') then
   g_cell%dis_posi(2)=gn%mdis_posi(1)
   select case(g_cell%edge(2))
    case(1,3); g_cell%dis_posi(2)=g_cell%dis_posi(2)+dfloat(gn%x_idx-g_cell%x_idx)
	case(2,4); g_cell%dis_posi(2)=g_cell%dis_posi(2)+dfloat(gn%y_idx-g_cell%y_idx)
   end select
  end if   	 
 end if
end if


end subroutine compute


subroutine normal_of_single(g_cell,normal)
! This subroutine computes the normal of a single critcal cell
! whose discontinuity curve segment is too short.

implicit none
type(geo_info), intent(in) :: g_cell
real(8), dimension(2), intent(out) :: normal

integer :: single

normal=error_data
call find_single(g_cell,single)

select case(single)
 case(1)
  normal(1)=0.5d0*dsqrt(2.0d0); normal(2)=0.5d0*dsqrt(2.0d0)
  if(g_cell%edge(1).eq.4) normal=-normal
 case(2)
  normal(1)=-0.5d0*dsqrt(2.0d0); normal(2)=0.5d0*dsqrt(2.0d0)
  if(g_cell%edge(1).eq.1) normal=-normal
 case(3)
  normal(1)=-0.5d0*dsqrt(2.0d0); normal(2)=-0.5d0*dsqrt(2.0d0)
  if(g_cell%edge(1).eq.2) normal=-normal
 case(4) 
  normal(1)=0.5d0*dsqrt(2.0d0); normal(2)=-0.5d0*dsqrt(2.0d0)
  if(g_cell%edge(1).eq.3) normal=-normal
end select

end subroutine normal_of_single


end subroutine compute_posi


subroutine compute_xy_mid_posi(acv)
! This subroutine recomputes the mdis_posi for $xy$-type discontinuity cells. The computation is 
! performed on a single discontinuity curve. The recomputed mdis_posi is the area ratio of the left 
! area computed from the mass of the ordinary state. The area ratio is stored in the first component
! of the mdis_posi. This information will be used in the numerical
! surface tension(regularization) on the discontinuity curve.

implicit none
type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
character*8 :: partner
integer :: head_mark, end_mark
logical :: head_switch

type(phy_info) :: p_cell

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
 g_cell=temp%g_cell
 if(g_cell%g_type.eq.'xy ') then
  p_cell=temp%p_cell
  g_cell%mdis_posi(1)=p_cell%or_state%value(1)-p_cell%r_state%value(1)
  g_cell%mdis_posi(1)=g_cell%mdis_posi(1)/(p_cell%l_state%value(1)-p_cell%r_state%value(1))
  g_cell%mdis_posi(2)=0.0d0
 end if
  
 temp%g_cell=g_cell

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine compute_xy_mid_posi


subroutine normal_of_corner(num,normal)
! This subroutine computes the normal for a corner point.

integer, intent(in) :: num
real(8), intent(out), dimension(2) :: normal

type(geo_info) :: g_cell, gn_cell

g_cell=temp%g_cell
select case(num)
 case(1); gn_cell=temp%p_nxt_dr%g_cell
 case(2); gn_cell=temp%n_nxt_dr%g_cell
 case default; call error_message
end select

if(iabs(g_cell%x_idx-gn_cell%x_idx)+iabs(g_cell%y_idx-gn_cell%y_idx).ne.1) then
 call error_message
end if

if(gn_cell%x_idx.gt.g_cell%x_idx) then
 normal=(/0.0d0,1.0d0/)
else
 if(gn_cell%x_idx.lt.g_cell%x_idx) then
  normal=(/0.0d0,-1.0d0/)
 else
  if(gn_cell%y_idx.gt.g_cell%y_idx) then
   normal=(/1.0d0,0.0d0/)
  else
   normal=(/-1.0d0,0.0d0/)
  end if
 end if
end if

if(num.eq.1) normal=-normal

end subroutine normal_of_corner


end module compute_discontinuity_positions