module manipulation_one

use solu_comput
! 'solu_com.f90'

use node_cells
! 'nd_cls.f90'

use variables_4_node_treat
! 'bs_var.f90'

implicit none

integer, parameter :: step_done=20

public  manipulation_1, manipulation_2, manipulation_3, &
        manipulation_5, manipulation_6, manipulation_8
private step_done, symmetric_mean_state, symmetric_mean_physical, symmetric_state, &
        symmetric_physical


contains


subroutine manipulation_1(nd_nmb)

implicit none

integer, intent(in) :: nd_nmb

type(state) :: or_state_ndd, difference
type(critical_cell), pointer :: temp

if(time_step.ge.2-step_done) then
 
 if(nd_nmb.eq.1) then
  call evaluate_node_state(1,or_state_ndd)
  difference=ndd(1)%n_cell%or_state-or_state_ndd
  ndd(1)%n_cell%or_state=or_state_ndd
  temp=>cvv(4)%eend%previous
  temp%p_cell%or_state=temp%p_cell%or_state+difference
 end if
 
 if(nd_nmb.eq.2) then
  call evaluate_node_state(2,or_state_ndd)
  difference=ndd(2)%n_cell%or_state-or_state_ndd
  ndd(2)%n_cell%or_state=or_state_ndd
  temp=>cvv(1)%eend%previous
  temp%p_cell%or_state=temp%p_cell%or_state+difference
 end if

end if

end subroutine manipulation_1


subroutine manipulation_2(recover_or_half_step)

implicit none

character*10, intent(in) :: recover_or_half_step

type(auxi_crit_cell), pointer :: temp, temp_23, temp_5, temp_14, &
                                 temp_7 
type(geo_info) :: g_cell, g_cell_5, g_cell_7, g_cell_14, g_cell_p, &
                  g_cell_77 
type(geo_info) :: ggg
real(8), dimension(2) :: position_1, position_2, pt_1, pt_2, &
                         difference_1, difference_2, difference, pt
integer :: x_1, y_1, x_2, y_2
real(8), dimension(2,2) :: pts1, pts2
character*10 :: err_mesg 
real(8) :: xx, rr

! The triple-point positions of the two node cells are computed as 
! the intersection points of the curve_1 and curve_5 and curve_4
! and curve_5
if(time_step.ge.2-step_done.and.time_step.le.10000-step_done) then
! First, compute the position in node_2.
 x_2=ndd(2)%n_cell%x_idx; y_2=ndd(2)%n_cell%y_idx
 temp_23=>acvv(2)%eend%previous
 pts1(1,1)=-1.0d0; pts1(1,2)=temp_23%g_cell%dis_posi(1)
 pts1(2,1)=1.0d0; pts1(2,2)=temp_23%g_cell%dis_posi(2)
 
 if(acvv(5)%begin_end.eq.2) then
  temp_5=>acvv(5)%begin%next
 else 
  if(acvv(5)%end_end.eq.2) then
   temp_5=>acvv(5)%eend%previous
  else
   call error_message 
  end if
 end if
    
 g_cell_5=temp_5%g_cell
 call pick_point(g_cell_5,1,pt_1); call pick_point(g_cell_5,2,pt_2)
 difference_1(1)=dfloat(x_2-g_cell_5%x_idx)
 difference_1(2)=dfloat(y_2-g_cell_5%y_idx)
 difference_2=difference_1
 rr=length_of_segment(pt_1,pt_2)
 if(rr.lt.0.1d0) then
  if(associated(temp_5%next)) then
   ggg=temp_5%next%g_cell
   call pick_point(ggg,2,pt_2)
   difference_2(1)=dfloat(x_2-ggg%x_idx)
   difference_2(2)=dfloat(y_2-ggg%y_idx)
  else
   pt_2=pt_1+(/1.0d0,-1.0d0/)
  end if
 end if
 pt_1=pt_1-difference_1; pt_2=pt_2-difference_2
 pts2(1,:)=pt_1; pts2(2,:)=pt_2
 
 call ints_of_lines(pts1,pts2,position_2,err_mesg)

! Second, compute the position in node_1.
 x_1=ndd(1)%n_cell%x_idx; y_1=ndd(1)%n_cell%y_idx
 temp_23=>acvv(3)%eend%previous
 pts1(1,1)=temp_23%g_cell%dis_posi(1); pts1(1,2)=-1.0d0
 pts1(2,1)=temp_23%g_cell%dis_posi(2); pts1(2,2)=1.0d0

 if(acvv(5)%begin_end.eq.1) then
  temp_5=>acvv(5)%begin%next
 else
  if(acvv(5)%end_end.eq.1) then
   temp_5=>acvv(5)%eend%previous
  else
   call error_message
  end if
 end if

 g_cell_5=temp_5%g_cell
 call pick_point(g_cell_5,1,pt_1); call pick_point(g_cell_5,2,pt_2)
 difference_2(1)=dfloat(x_1-g_cell_5%x_idx)
 difference_2(2)=dfloat(y_1-g_cell_5%y_idx)
 difference_1=difference_2
 rr=length_of_segment(pt_1,pt_2)
 if(rr.lt.0.1d0) then
  if(associated(temp_5%previous)) then
   ggg=temp_5%previous%g_cell
   call pick_point(ggg,1,pt_1)
   difference_1(1)=dfloat(x_1-ggg%x_idx)
   difference_1(2)=dfloat(y_1-ggg%y_idx)
  else
   pt_2=pt_1+(/1.0d0,-1.0d0/)
  end if    
 end if
 pt_1=pt_1-difference_1; pt_2=pt_2-difference_2
 pts2(1,:)=pt_1; pts2(2,:)=pt_2

 call ints_of_lines(pts1,pts2,position_1,err_mesg)

 if(time_step.le.1000-step_done) then
  xx=0.5d0*(position_2(1)+position_1(2))
  position_2(1)=xx; position_1(2)=xx
 end if
 
 temp_14=>acvv(1)%eend%previous

 g_cell_14=temp_14%g_cell
 if(g_cell_14%g_type.eq.'xx ') then
  select case(recover_or_half_step)
   case('recover   '); pts1(1,1)=g_cell_14%mdis_posi(1); pts1(1,2)=0.0d0
   case('half_step '); call pick_point(g_cell_14,1,pt_1); pts1(1,:)=pt_1
   case default; call error_message
  end select
 else
  if(g_cell_14%g_type.eq.'xy ') then
   select case(recover_or_half_step)
    case('recover   '); call pick_middle_point(g_cell_14,1,pt_1)
    case('half_step '); call pick_point(g_cell_14,1,pt_1)
    case default; call error_message
   end select
   pts1(1,:)=pt_1
  else
   select case(recover_or_half_step)
    case('recover   '); pts1(1,1)=0.0d0; pts1(1,2)=g_cell_14%mdis_posi(2)
	case('half_step '); call pick_point(g_cell_14,1,pt_1); pts1(1,:)=pt_1
	case default; call error_message
   end select
  end if
 end if
 difference=(/dfloat(x_2-g_cell_14%x_idx),dfloat(y_2-g_cell_14%y_idx)/)
 pts1(2,:)=position_2+difference
 select case(g_cell_14%edge(2))
  case(1); pts2(1,:)=-0.5d0; pts2(2,:)=(/0.5d0,-0.5d0/)
  case(4); pts2(1,:)=-0.5d0; pts2(2,:)=(/-0.5d0,0.5d0/)
  case default; call error_message
 end select
 call ints_of_lines(pts1,pts2,pt,err_mesg)
 select case(g_cell_14%edge(2))
  case(1); g_cell_14%dis_posi(2)=pt(1)
  case(4); g_cell_14%dis_posi(2)=pt(2)
 end select
 temp_14%g_cell=g_cell_14
 if(temp_14%partner.eq.'next    ') call error_message
 if(temp_14%partner.eq.'previous') then
  g_cell_p=temp_14%previous%g_cell
  difference=(/dfloat(g_cell_p%x_idx-g_cell_14%x_idx),dfloat(g_cell_p%y_idx-g_cell_14%y_idx)/)
  temp_14%previous%g_cell%dis_posi(2)=g_cell_14%dis_posi(2)-difference(1)-difference(2)
 end if

 temp_14=>acvv(4)%eend%previous
 g_cell_14=temp_14%g_cell
 if(g_cell_14%g_type.eq.'xx ') then
  select case('recover_or_half_step')
   case('recover   '); pts1(1,1)=g_cell_14%mdis_posi(1); pts1(1,2)=0.0d0
   case('half_step '); call pick_point(g_cell_14,1,pt_1); pts1(1,:)=pt_1
   case default; call error_message
  end select
 else
  if(g_cell_14%g_type.eq.'xy ') then
   select case(recover_or_half_step) 
    case('recover   '); call pick_middle_point(g_cell_14,1,pt_1)
    case('half_step '); call pick_point(g_cell_14,1,pt_1)
	case default; call error_message
   end select	 
   pts1(1,:)=pt_1
  else
   select case(recover_or_half_step)
    case('recover   '); pts1(1,1)=0.0d0; pts1(1,2)=g_cell_14%mdis_posi(2)
    case('half_step '); call pick_point(g_cell_14,1,pt_1); pts1(1,:)=pt_1
    case default; call error_message
   end select
  end if
 end if
 difference=(/dfloat(x_1-g_cell_14%x_idx),dfloat(y_1-g_cell_14%y_idx)/)
 pts1(2,:)=position_1+difference
 select case(g_cell_14%edge(2))
  case(1); pts2(1,:)=-0.5d0; pts2(2,:)=(/0.5d0,-0.5d0/)
  case(4); pts2(1,:)=-0.5d0; pts2(2,:)=(/-0.5d0,0.5d0/)
  case default; call error_message
 end select
 call ints_of_lines(pts1,pts2,pt,err_mesg)
 select case(g_cell_14%edge(2))
  case(1); g_cell_14%dis_posi(2)=pt(1)
  case(4); g_cell_14%dis_posi(2)=pt(2)
 end select
 temp_14%g_cell=g_cell_14
 if(temp_14%partner.eq.'next    ') call error_message
 if(temp_14%partner.eq.'previous') then
  g_cell_p=temp_14%previous%g_cell
  difference=(/dfloat(g_cell_p%x_idx-g_cell_14%x_idx),dfloat(g_cell_p%y_idx-g_cell_14%y_idx)/)
  temp_14%previous%g_cell%dis_posi(2)=g_cell_14%dis_posi(2)-difference(1)-difference(2)
 end if
 
end if

end subroutine manipulation_2


subroutine manipulation_3

implicit none

type(auxi_crit_cell), pointer :: temp

if(time_step.eq.15-step_done) then
 temp=>acvv(5)%eend%previous
 temp%g_cell%dis_posi=-0.5000001d0
 temp=>temp%previous
 temp%g_cell%dis_posi(2)=-0.5000001d0
 temp=>temp%previous
 temp%g_cell%dis_posi(2)=-1.5000001d0
end if

if(time_step.eq.34-step_done) then
 temp=>acvv(5)%eend%previous
 temp%g_cell%dis_posi=-0.5000001d0
 temp=>temp%previous
 temp%g_cell%dis_posi(2)=-0.5000001d0
 temp=>temp%previous
 temp%g_cell%dis_posi(2)=-1.5000001d0
end if

if(time_step.eq.44-step_done) then
 temp=>acvv(5)%eend%previous
 temp%g_cell%dis_posi(2)=-0.5000001d0
 temp=>temp%previous
 temp%g_cell%dis_posi(2)=-1.5000001d0
end if 

if(time_step.eq.55-step_done) then
 temp=>acvv(5)%begin%next
 temp%g_cell%dis_posi=-0.5000001d0
 temp=>temp%next
 temp%g_cell%dis_posi(1)=-0.5000001d0
 temp=>temp%next
 temp%g_cell%dis_posi(1)=-1.5000001d0
end if

if(time_step.eq.66-step_done) then
 temp=>acvv(5)%eend%previous
 temp%g_cell%dis_posi(2)=-0.5000001d0
 temp=>temp%previous
 temp%g_cell%dis_posi(2)=-1.5000001d0
end if

if(time_step.eq.93-step_done) then
 temp=>acvv(5)%begin%next
 temp%g_cell%dis_posi=-0.5000001d0
 temp=>temp%next
 temp%g_cell%dis_posi(1)=-0.5000001d0
 temp=>temp%next
 temp%g_cell%dis_posi(1)=-1.5000001d0
end if

end subroutine manipulation_3


subroutine manipulation_31

implicit none

type(state) :: state_1, state_2

if(time_step.ge.15-step_done) then
 call symmetric_mean_state(uu(-1,0),uu(0,-1),state_1,state_2)
 uu(-1,0)=state_1
 uu(0,-1)=state_2
end if

if(time_step.ge.34-step_done) then
 call symmetric_mean_state(uu(-2,-1),uu(-1,-2),state_1,state_2)
 uu(-2,-1)=state_1
 uu(-1,-2)=state_2
end if

if(time_step.ge.55-step_done) then
 call symmetric_mean_state(uu(-3,-2),uu(-2,-3),state_1,state_2)
 uu(-3,-2)=state_1
 uu(-2,-3)=state_2
end if

if(time_step.ge.76-step_done) then
 call symmetric_mean_state(uu(-4,-3),uu(-3,-4),state_1,state_2)
 uu(-4,-3)=state_1
 uu(-3,-4)=state_2
end if

if(time_step.ge.85-step_done) then
 call symmetric_mean_state(uu(-5,-3),uu(-3,-5),state_1,state_2)
 uu(-5,-3)=state_1
 uu(-3,-5)=state_2
end if

if(time_step.ge.93-step_done) then
 call symmetric_mean_state(uu(-5,-4),uu(-4,-5),state_1,state_2)
 uu(-5,-4)=state_1
 uu(-4,-5)=state_2
end if

end subroutine manipulation_31


subroutine manipulation_4

implicit none

type(critical_cell), pointer :: temp
real(8) :: area_l, area_r
type(state) :: state_new, difference

temp=>cvv(5)%begin%next
area_l=side_area(temp%g_cell,'left  ')
area_r=side_area(temp%g_cell,'right ')
state_new=area_l*temp%p_cell%l_state+area_r*temp%p_cell%r_state
difference=temp%p_cell%or_state-state_new
temp%p_cell%or_state=state_new
ndd(cvv(5)%begin_end)%n_cell%or_state=ndd(cvv(5)%begin_end)%n_cell%or_state+difference

temp=>cvv(5)%eend%previous
area_l=side_area(temp%g_cell,'left  ')
area_r=side_area(temp%g_cell,'right ')
state_new=area_l*temp%p_cell%l_state+area_r*temp%p_cell%r_state
difference=temp%p_cell%or_state-state_new
temp%p_cell%or_state=state_new
ndd(cvv(5)%end_end)%n_cell%or_state=ndd(cvv(5)%end_end)%n_cell%or_state+difference

end subroutine manipulation_4


subroutine manipulation_5

implicit none

type(curve_plug), pointer :: tempp_1, tempp_2
type(cv_plug_info) :: plug_1, plug_2
type(state), dimension(-1:1,-1:1) :: uu_1, uu_2, uuu1, uuu2
integer :: i, j
logical :: if_do
character*6 :: node1_5, node2_5

if_do=.false.

!if(time_step.eq.-step_done) if_do=.true.
if(time_step.ge.2-step_done.and.time_step.le.1000-step_done) if_do=.true.

if(if_do) then
 select case(acvv(5)%begin_end)
  case(1); node1_5='begin '
  case(2); node2_5='begin '
  case default; call error_message
 end select
 select case(acvv(5)%end_end)
  case(1); node1_5='end   '
  case(2); node2_5='end   '
  case default; call error_message
 end select
 if(node1_5.eq.node2_5) call error_message

 call visit(1,5,node1_5,tempp_1); plug_1=tempp_1%plug
 call visit(2,5,node2_5,tempp_2); plug_2=tempp_2%plug

 do i=-1,1; do j=-1,1
   uu_1(i,j)=plug_1%u_front(i,j)
   uu_2(i,j)=plug_2%u_behind(i,j)
  end do; end do

 do i=-1,1; do j=-1,1
  uuu1(i,j)%value(1)=0.5d0*(uu_1(i,j)%value(1)+uu_2(j,i)%value(1))
  uuu1(i,j)%value(2)=0.5d0*(uu_1(i,j)%value(2)+uu_2(j,i)%value(3))
  uuu1(i,j)%value(3)=0.5d0*(uu_1(i,j)%value(3)+uu_2(j,i)%value(2))
  uuu1(i,j)%value(4)=0.5d0*(uu_1(i,j)%value(4)+uu_2(j,i)%value(4))
  uuu2(i,j)%value(1)=0.5d0*(uu_1(i,j)%value(1)+uu_2(j,i)%value(1))
  uuu2(i,j)%value(2)=0.5d0*(uu_1(i,j)%value(3)+uu_2(j,i)%value(2))
  uuu2(i,j)%value(3)=0.5d0*(uu_1(i,j)%value(2)+uu_2(j,i)%value(3))
  uuu2(i,j)%value(4)=0.5d0*(uu_1(i,j)%value(4)+uu_2(j,i)%value(4))              
 end do; end do

 do i=-1,1; do j=-1,1
  plug_1%u_front(i,j)=uuu1(i,j)
  plug_2%u_behind(i,j)=uuu2(j,i)
 end do; end do

 tempp_1%plug=plug_1; tempp_2%plug=plug_2
  
 do i=-1,1; do j=-1,1
  tempp_1%next%plug%u_behind(i,j)=plug_1%u_front(i,j)
  tempp_2%previous%plug%u_front(i,j)=plug_2%u_behind(i,j)
 end do; end do

 do i=-1,1; do j=-1,1
   uu_1(i,j)=plug_1%u_behind(i,j)
   uu_2(i,j)=plug_2%u_front(i,j)
  end do; end do

 do i=-1,1; do j=-1,1
  uuu1(i,j)%value(1)=0.5d0*(uu_1(i,j)%value(1)+uu_2(j,i)%value(1))
  uuu1(i,j)%value(2)=0.5d0*(uu_1(i,j)%value(2)+uu_2(j,i)%value(3))
  uuu1(i,j)%value(3)=0.5d0*(uu_1(i,j)%value(3)+uu_2(j,i)%value(2))
  uuu1(i,j)%value(4)=0.5d0*(uu_1(i,j)%value(4)+uu_2(j,i)%value(4))
  uuu2(i,j)%value(1)=0.5d0*(uu_1(i,j)%value(1)+uu_2(j,i)%value(1))
  uuu2(i,j)%value(2)=0.5d0*(uu_1(i,j)%value(3)+uu_2(j,i)%value(2))
  uuu2(i,j)%value(3)=0.5d0*(uu_1(i,j)%value(2)+uu_2(j,i)%value(3))
  uuu2(i,j)%value(4)=0.5d0*(uu_1(i,j)%value(4)+uu_2(j,i)%value(4))              
 end do; end do

 do i=-1,1; do j=-1,1
  plug_1%u_behind(i,j)=uuu1(i,j)
  plug_2%u_front(i,j)=uuu2(j,i)
 end do; end do
  
 tempp_1%plug=plug_1; tempp_2%plug=plug_2
  
 do i=-1,1; do j=-1,1
  tempp_1%previous%plug%u_front(i,j)=plug_1%u_behind(i,j)
  tempp_2%next%plug%u_behind(i,j)=plug_2%u_front(i,j)
 end do; end do

end if

end subroutine manipulation_5


subroutine manipulation_6

implicit none

type(state) :: uu_1, uu_2
logical :: if_do

if_do=.false.

!if(time_step.eq.34-step_done)  if_do=.true.
if(time_step.ge.2-step_done.and.time_step.le.1000-step_done) if_do=.true.

if(if_do) then
  uu_1%value(1)=0.5d0*(ndd(1)%n_cell%or_state%value(1)+ndd(2)%n_cell%or_state%value(1))
  uu_1%value(2)=0.5d0*(ndd(1)%n_cell%or_state%value(2)+ndd(2)%n_cell%or_state%value(3))
  uu_1%value(3)=0.5d0*(ndd(1)%n_cell%or_state%value(3)+ndd(2)%n_cell%or_state%value(2))
  uu_1%value(4)=0.5d0*(ndd(1)%n_cell%or_state%value(4)+ndd(2)%n_cell%or_state%value(4))
  uu_2%value(1)=0.5d0*(ndd(1)%n_cell%or_state%value(1)+ndd(2)%n_cell%or_state%value(1))
  uu_2%value(2)=0.5d0*(ndd(1)%n_cell%or_state%value(3)+ndd(2)%n_cell%or_state%value(2))
  uu_2%value(3)=0.5d0*(ndd(1)%n_cell%or_state%value(2)+ndd(2)%n_cell%or_state%value(3))
  uu_2%value(4)=0.5d0*(ndd(1)%n_cell%or_state%value(4)+ndd(2)%n_cell%or_state%value(4))
  ndd(1)%n_cell%or_state=uu_1; ndd(2)%n_cell%or_state=uu_2
end if

end subroutine manipulation_6


subroutine manipulation_7

implicit none

type(critical_cell), pointer :: temp
real(8) :: area_l, area_r
type(state) :: state_new, difference
logical :: if_do

if(cvv(5)%status.ne.'awake ') return

if(cvv(5)%total.le.1) then
 
 if_do=.false.
 temp=>cvv(5)%begin%next
 if(cvv(5)%begin_end.eq.2.and.temp%g_cell%edge(1).eq.4) if_do=.true.
 if(cvv(5)%end_end.eq.2.and.temp%g_cell%edge(2).eq.4) if_do=.true.
 if(if_do) then
  area_l=side_area(temp%g_cell,'left  ')
  area_r=side_area(temp%g_cell,'right ')
  state_new=area_l*temp%p_cell%l_state+area_r*temp%p_cell%r_state
  difference=temp%p_cell%or_state-state_new
  temp%p_cell%or_state=state_new
  ndd(cvv(5)%begin_end)%n_cell%or_state=ndd(cvv(5)%begin_end)%n_cell%or_state+0.5d0*difference
  ndd(cvv(5)%end_end)%n_cell%or_state=ndd(cvv(5)%end_end)%n_cell%or_state+0.5d0*difference
 end if

else

 if_do=.false.
 temp=>cvv(5)%begin%next
 if(cvv(5)%begin_end.eq.1.and.temp%g_cell%edge(1).eq.1) if_do=.true.
 if(cvv(5)%begin_end.eq.2.and.temp%g_cell%edge(1).eq.4) if_do=.true. 
 if(if_do) then
  area_l=side_area(temp%g_cell,'left  ')
  area_r=side_area(temp%g_cell,'right ')
  state_new=area_l*temp%p_cell%l_state+area_r*temp%p_cell%r_state
  difference=temp%p_cell%or_state-state_new
  temp%p_cell%or_state=state_new
  ndd(cvv(5)%begin_end)%n_cell%or_state=ndd(cvv(5)%begin_end)%n_cell%or_state+difference
 end if
 
 if_do=.false. 
 temp=>cvv(5)%eend%previous
 if(cvv(5)%end_end.eq.1.and.temp%g_cell%edge(2).eq.1) if_do=.true.
 if(cvv(5)%end_end.eq.2.and.temp%g_cell%edge(2).eq.4) if_do=.true. 
 if(if_do) then
  area_l=side_area(temp%g_cell,'left  ')
  area_r=side_area(temp%g_cell,'right ')
  state_new=area_l*temp%p_cell%l_state+area_r*temp%p_cell%r_state
  difference=temp%p_cell%or_state-state_new
  temp%p_cell%or_state=state_new
  ndd(cvv(5)%end_end)%n_cell%or_state=ndd(cvv(5)%end_end)%n_cell%or_state+difference

 end if

end if

end subroutine manipulation_7


subroutine manipulation_8
! This manipulation is to maintain the symmetry of the computation
! at the begining stage.

implicit none

integer :: ii
type(auxi_crit_cell), pointer :: temp_b, temp_e, temp_g
type(geo_info) :: gc_b, gc_e
type(phy_info) :: pc_b, pc_bb, pc_e, pc_ee
type(state) :: state_1, state_2, difference
real(8) :: posi_1, posi_2, area_l, area_r

if(time_step.ge.2-step_done.and.time_step.le.1000-step_done) then

! The symmetry of curve_5.
 temp_b=>acvv(5)%begin%next; temp_e=>acvv(5)%eend%previous
 gc_b=temp_b%g_cell; gc_e=temp_e%g_cell
 do while(gc_b%x_idx.ne.gc_b%y_idx)
  if(gc_b%x_idx.ne.gc_e%y_idx.or.gc_b%y_idx.ne.gc_e%x_idx) call error_message
  pc_b=temp_b%p_cell; pc_e=temp_e%p_cell
  posi_1=0.5d0*(gc_b%dis_posi(1)+gc_e%dis_posi(2))
  posi_2=0.5d0*(gc_b%dis_posi(2)+gc_e%dis_posi(1))
  gc_b%dis_posi(1)=posi_1; gc_b%dis_posi(2)=posi_2
  gc_e%dis_posi(1)=posi_2; gc_e%dis_posi(2)=posi_1
  call symmetric_mean_physical(pc_b,pc_e,pc_bb,pc_ee)
  temp_b%g_cell=gc_b; temp_e%g_cell=gc_e
  temp_b%p_cell=pc_bb; temp_e%p_cell=pc_ee
  temp_b=>temp_b%next; temp_e=>temp_e%previous
  gc_b=temp_b%g_cell; gc_e=temp_e%g_cell
 end do
 if(gc_b%x_idx.ne.gc_e%y_idx.or.gc_b%y_idx.ne.gc_e%x_idx) call error_message
 posi_1=0.5d0*(gc_b%dis_posi(1)+gc_e%dis_posi(2))
 gc_b%dis_posi(1)=posi_1; gc_b%dis_posi(2)=posi_1
 pc_b=temp_b%p_cell
 call symmetric_physical(pc_b,pc_bb)
 temp_b%g_cell=gc_b; temp_b%p_cell=pc_bb

! The symmetry of curve_1 and curve_4
 temp_b=>acvv(1)%eend%previous; temp_e=>acvv(4)%eend%previous
 do ii=1,2
  gc_b=temp_b%g_cell; gc_e=temp_e%g_cell
  pc_b=temp_b%p_cell; pc_e=temp_e%p_cell
  posi_1=0.5d0*(gc_b%dis_posi(1)+gc_e%dis_posi(1))
  posi_2=0.5d0*(gc_b%dis_posi(2)+gc_e%dis_posi(2))
  gc_b%dis_posi(1)=posi_1; gc_b%dis_posi(2)=posi_2
  gc_e%dis_posi(1)=posi_1; gc_e%dis_posi(2)=posi_2
  call symmetric_mean_state(pc_b%l_state,pc_e%r_state,state_1,state_2)
  pc_b%l_state=state_1; pc_e%r_state=state_2
  call symmetric_mean_state(pc_b%r_state,pc_e%l_state,state_1,state_2)
  pc_b%r_state=state_1; pc_e%l_state=state_2
  call symmetric_mean_state(pc_b%or_state,pc_e%or_state,state_1,state_2)
  pc_b%or_state=state_1; pc_e%or_state=state_2
  temp_b%g_cell=gc_b; temp_b%p_cell%or_state=state_1
  temp_e%g_cell=gc_e; temp_e%p_cell%or_state=state_2

  temp_b=>temp_b%previous; temp_e=>temp_e%previous

 end do

end if

end subroutine manipulation_8


subroutine symmetric_mean_state(state_1,state_2,state_m1,state_m2)

implicit none

type(state), intent(in) :: state_1, state_2
type(state), intent(out) :: state_m1, state_m2

state_m1%value(1)=0.5d0*(state_1%value(1)+state_2%value(1))
state_m1%value(2)=0.5d0*(state_1%value(2)+state_2%value(3))
state_m1%value(3)=0.5d0*(state_1%value(3)+state_2%value(2))
state_m1%value(4)=0.5d0*(state_1%value(4)+state_2%value(4))
state_m2%value(1)=0.5d0*(state_1%value(1)+state_2%value(1))
state_m2%value(3)=0.5d0*(state_1%value(2)+state_2%value(3))
state_m2%value(2)=0.5d0*(state_1%value(3)+state_2%value(2))
state_m2%value(4)=0.5d0*(state_1%value(4)+state_2%value(4))

end subroutine symmetric_mean_state


subroutine symmetric_mean_physical(p_cell1,p_cell2,p_cellm1,p_cellm2)

implicit none

type(phy_info), intent(in) :: p_cell1, p_cell2
type(phy_info), intent(out) :: p_cellm1, p_cellm2

p_cellm1%wv_nb=p_cell1%wv_nb; p_cellm2%wv_nb=p_cell2%wv_nb

call symmetric_mean_state(p_cell1%l_state,p_cell2%l_state,p_cellm1%l_state,p_cellm2%l_state)
call symmetric_mean_state(p_cell1%r_state,p_cell2%r_state,p_cellm1%r_state,p_cellm2%r_state)
call symmetric_mean_state(p_cell1%or_state,p_cell2%or_state,p_cellm1%or_state,p_cellm2%or_state)

end subroutine symmetric_mean_physical


subroutine symmetric_state(state_i,state_m)

implicit none

type(state), intent(in) :: state_i
type(state), intent(out) :: state_m

state_m=state_i
state_m%value(2)=0.5d0*(state_i%value(2)+state_i%value(3))
state_m%value(3)=0.5d0*(state_i%value(3)+state_i%value(2))

end subroutine symmetric_state


subroutine symmetric_physical(p_cell,p_cellm)

implicit none

type(phy_info), intent(in) :: p_cell
type(phy_info), intent(out) :: p_cellm

p_cellm%wv_nb=p_cell%wv_nb

call symmetric_state(p_cell%l_state,p_cellm%l_state)
call symmetric_state(p_cell%r_state,p_cellm%r_state)
call symmetric_state(p_cell%or_state,p_cellm%or_state)

end subroutine symmetric_physical


subroutine symmetric_mean_geometry(g_cell1,g_cell2,g_cellm1,g_cellm2)

implicit none

type(geo_info), intent(in) :: g_cell1, g_cell2
type(geo_info), intent(out) :: g_cellm1, g_cellm2

real(8) :: xx

g_cellm1=g_cell1; g_cellm2=g_cell2

xx=0.5d0*(g_cell1%dis_posi(1)+g_cell2%dis_posi(2))
g_cellm1%dis_posi(1)=xx; g_cellm2%dis_posi(2)=xx
xx=0.5d0*(g_cell1%dis_posi(2)+g_cell2%dis_posi(1))
g_cellm1%dis_posi(2)=xx; g_cellm2%dis_posi(1)=xx

end subroutine symmetric_mean_geometry


subroutine symmetric_geometry(g_cell,g_cellm)

implicit none

type(geo_info), intent(in) :: g_cell
type(geo_info), intent(out) :: g_cellm

real(8) :: xx

g_cellm=g_cell

xx=0.5d0*(g_cell%dis_posi(1)+g_cell%dis_posi(2))
g_cellm%dis_posi(1)=xx; g_cellm%dis_posi(2)=xx

end subroutine symmetric_geometry


end module manipulation_one