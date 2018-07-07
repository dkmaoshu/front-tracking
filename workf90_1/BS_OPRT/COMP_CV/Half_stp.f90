module half_step_computation
! This subroutine is for computing temporary discontinuity 
! positions, which will be used later in evaluating of numerical
! fluxes, using Euler forward method.  Actually, we compute the
! discontinuity positions at the next time and store the information
! in 'g_cell%dis_posi'. The current discontinuity positions are
! then stored in 'c_cell%tdp'. The half-step discontinuity positions
! are then the mean values of the two positions.

use auxiliary_discontinuity_curves
!'auxi_cv.f90'

implicit none


public  half_step
private half_step_operation



contains


subroutine half_step

implicit none

integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake') then
  call store_old_positions(acvv(i))
  call half_step_operation(acvv(i))
! Carry out the half-step computation.

!  call check_list_ac(1,'down  ')

 end if
end do

end subroutine half_step


subroutine store_old_positions(acv)

implicit none
type(auxi_discv), intent(in) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(adss_info) :: address
type(geo_info) :: g_cell
type(comp_info) :: c_cell
integer :: head_mark, end_mark
logical :: head_switch

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

 g_cell=temp%g_cell; c_cell=temp%c_cell
 address=temp%address

! Store the original discontinuity position in 'c_cell%tdp'.
 c_cell%tdp=g_cell%dis_posi

 temp%g_cell=g_cell; temp%c_cell=c_cell 
  
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine store_old_positions


subroutine half_step_operation(acv)
! Carry out the half-step computation in an individual 
! discontinuity curve.

implicit none
type(auxi_discv), intent(in) :: acv

type(auxi_crit_cell), pointer :: temp, boundary_end
type(adss_info) :: address
type(geo_info) :: g_cell
type(comp_info) :: c_cell
character*8 :: partner
integer :: head_mark, end_mark
logical :: head_switch

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

 g_cell=temp%g_cell
 partner=temp%partner; address=temp%address

 call computation

 temp%g_cell=g_cell 
  
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do


contains


subroutine computation

implicit none

type(auxi_crit_cell), pointer :: tempap
type(geo_info) :: gp_cell
integer :: i

select case(partner)

 case('previous')
! Transfer the discontinuity positions from the previous partner.
  tempap=>temp%previous; gp_cell=tempap%g_cell
  g_cell%dis_posi=gp_cell%dis_posi
  do i=1,2
   select case(g_cell%edge(i))
    case(1,3)
     g_cell%dis_posi(i)=g_cell%dis_posi(i)+ &
	                    dfloat(gp_cell%x_idx-g_cell%x_idx)
    case(2,4)
     g_cell%dis_posi(i)=g_cell%dis_posi(i)+ &
	                    dfloat(gp_cell%y_idx-g_cell%y_idx)
   end select
  end do
  temp%g_cell=g_cell
 
 case('next    ','single  ')
  if(.not.associated(temp%p_nxt_dr)) call compute(1)
  if(associated(temp)) then
   call compute(2)
   temp%g_cell=g_cell
   if(associated(temp%n_nxt_dr)) then
    tempap=>temp%n_nxt_dr; gp_cell=tempap%g_cell
    gp_cell%dis_posi(1)=g_cell%dis_posi(2)
    gp_cell%normal(1,:)=g_cell%normal(2,:)
    select case(gp_cell%edge(1))
     case(1,3)
      gp_cell%dis_posi(1)=gp_cell%dis_posi(1)+ &
	                      dfloat(g_cell%x_idx-gp_cell%x_idx)
     case(2,4)
      gp_cell%dis_posi(1)=gp_cell%dis_posi(1)+ &
	                      dfloat(g_cell%y_idx-gp_cell%y_idx)
    end select
    tempap%g_cell=gp_cell
   end if
  end if
end select


end subroutine computation


subroutine compute(num)
! Compute the new discontinuity position using Euler forward method.

implicit none
integer, intent(in) :: num

type(geo_info) :: gnb_cell
type(phy_info) :: p_cell
type(comp_info) ::cnb_cell
type(state) :: uu_l, uu_r, uul, uur
real(8), dimension(2) :: normal, ppt
logical :: if_neighbored
character*12 :: wave_type
real(8) :: speed, displacement, x_posi, y_posi

call clean_up_g_cell(gnb_cell); call clean_up_c_cell(cnb_cell)
call clean_up_c_cell(c_cell)
if_neighbored=.true.

select case(num)

 case(1)
  c_cell=temp%c_cell
  if(associated(temp%p_nxt_dr)) then
   gnb_cell=temp%p_nxt_dr%g_cell; cnb_cell=temp%p_nxt_dr%c_cell
  else
   if_neighbored=.false.
  end if
  call pick_point(g_cell,1,ppt)

 case(2)
  select case(partner)
   case('single  '); c_cell=temp%c_cell
   case('next    '); c_cell=temp%next%c_cell
   case default; call error_message
  end select
  if(associated(temp%n_nxt_dr)) then
   gnb_cell=temp%n_nxt_dr%g_cell; cnb_cell=temp%n_nxt_dr%c_cell
  else
   if_neighbored=.false.
  end if
  call pick_point(g_cell,2,ppt)

end select

if(if_neighbored) then
 uu_l=0.5d0*(c_cell%l_state(0)+cnb_cell%l_state(0))
 uu_r=0.5d0*(c_cell%r_state(0)+cnb_cell%r_state(0))
else
 uu_l=c_cell%l_state(0); uu_r=c_cell%r_state(0)
end if

normal=g_cell%normal(num,:); p_cell=temp%p_cell
x_posi=(ppt(1)+dfloat(g_cell%x_idx))*h
y_posi=(ppt(2)+dfloat(g_cell%y_idx))*h
call riemann(uu_l,uu_r,normal,p_cell%wv_nb,x_posi,y_posi,uul,uur,speed,wave_type)

select case(g_cell%edge(num))
 case(1,3); speed=speed/normal(1)
 case(2,4); speed=speed/normal(2)
end select

displacement=speed*r
g_cell%dis_posi(num)=g_cell%dis_posi(num)+displacement


end subroutine compute

end subroutine half_step_operation


end module half_step_computation