module update_left_n_right_states
! This module is for updating the left and right states in critical cells after the 
! computation of ordinary states.

use evaluate_updated_xxyy
! 'eva_sts.f90'

implicit none

public  update_left_n_right
private update_xx_yy, update_xy, g_cell, p_cell, c_cell


contains


subroutine update_left_n_right

integer :: i

do i=1,cvn
 if(acvv(i)%status.eq.'awake ') call update_xx_yy(acvv(i))
end do

! call check_list_ac(1,'down  ')

do i=1,cvn
 if(acvv(i)%status.eq.'awake ') call return_data(acvv(i),1)
end do

! call check_partner_conserv_ac(1)

! call check_list_ac(1,'down  ')

do i=1,cvn
 if(acvv(i)%status.eq.'awake ') call update_xy(acvv(i))
end do

! call check_partner_conserv_ac(1)

do i=1,cvn
 if(acvv(i)%status.eq.'awake ') call return_data(acvv(i),2)
end do    

! call check_partner_conserv_ac(1)

end subroutine update_left_n_right


subroutine return_data(acv,no)
! Return the data stored in 'c_cell' back to the left and right states in 'xx'- and 'yy'-type
! critical cells.
implicit none
type(auxi_discv), intent(inout) :: acv
integer, intent(in) :: no

type(adss_info) :: address !, address1
integer :: head_mark, end_mark
logical :: head_switch
type(state) :: difference

type(phy_info) :: fff
type(state) :: ggg

real(8) :: xx, yy, rr
real(8) :: u_radial_l, u_radial_r, u_tangen_l, u_tangen_r

type(geo_info) :: gptn_cell
type(phy_info) :: pptn_cell

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
 g_cell=temp%g_cell; p_cell=temp%p_cell; c_cell=temp%c_cell
 address=temp%address
 fff=p_cell
 ggg=c_cell%l_state(0)
 select case(no)
  case(1)
   if(g_cell%g_type.eq.'xx '.or.g_cell%g_type.eq.'yy ') then
    difference=c_cell%l_state(0)-p_cell%l_state
    p_cell%l_state=c_cell%l_state(0)
    
    c_cell%l_state(0)=difference
    
	call update_stacked_neighbor('left  ')
	difference=c_cell%r_state(0)-p_cell%r_state
    p_cell%r_state=c_cell%r_state(0)
    
    c_cell%r_state(0)=difference
    
	ggg=p_cell%r_state
	call update_stacked_neighbor('right ')
   end if
  case(2)
   if(g_cell%g_type.eq.'xy ') then
    difference=c_cell%l_state(0)-p_cell%l_state
    p_cell%l_state=c_cell%l_state(0)
    
    c_cell%l_state(0)=difference
    
	call update_stacked_neighbor('left  ')
	difference=c_cell%r_state(0)-p_cell%r_state
    p_cell%r_state=c_cell%r_state(0)
    
    c_cell%r_state(0)=difference
    
	call update_stacked_neighbor('right ')
   end if
 end select

 temp%p_cell=p_cell; temp%c_cell=c_cell

 temp=>temp%next !; address1=temp%address
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

if(no.eq.1) then

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
  g_cell=temp%g_cell; p_cell=temp%p_cell; c_cell=temp%c_cell
  
  if(g_cell%g_type.eq.'xx '.or.g_cell%g_type.eq.'yy ') then
   
   select case(temp%partner)
    
    case('previous')
     gptn_cell=temp%previous%g_cell; pptn_cell=temp%previous%p_cell
     if(g_cell%x_idx+g_cell%y_idx.gt.gptn_cell%x_idx+gptn_cell%y_idx) then
      if(gptn_cell%point(3).eq.'left  ') then
       pptn_cell%or_state=pptn_cell%or_state-c_cell%l_state(0)
      else
       pptn_cell%or_state=pptn_cell%or_state-c_cell%r_state(0)
      end if	  
     else
      if(gptn_cell%point(1).eq.'left  ') then
       pptn_cell%or_state=pptn_cell%or_state-c_cell%l_state(0)
      else
       pptn_cell%or_state=pptn_cell%or_state-c_cell%r_state(0)
      end if	 
     end if
     temp%previous%p_cell=pptn_cell 
      	 	 	  	    
    case('next    ')
     gptn_cell=temp%next%g_cell; pptn_cell=temp%next%p_cell
     if(g_cell%x_idx+g_cell%y_idx.gt.gptn_cell%x_idx+gptn_cell%y_idx) then
      if(gptn_cell%point(3).eq.'left  ') then
       pptn_cell%or_state=pptn_cell%or_state-c_cell%l_state(0)
      else
       pptn_cell%or_state=pptn_cell%or_state-c_cell%r_state(0)
      end if	  
     else
      if(gptn_cell%point(1).eq.'left  ') then
       pptn_cell%or_state=pptn_cell%or_state-c_cell%l_state(0)
      else
       pptn_cell%or_state=pptn_cell%or_state-c_cell%r_state(0)
      end if	 
     end if
     temp%next%p_cell=pptn_cell

    case('single  ')
     
    case default; call error_message
    
   end select
  
  end if
   
  temp=>temp%next !; address1=temp%address
  if(associated(temp)) then
   head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
  else
   exit
  end if
 end do

end if


contains


subroutine update_stacked_neighbor(side)
! This subroutine updates the left or right states in neighboring critical cells which are
! stacked in the same grid cell. It also updates the ordinary cell-average in the present 
! critical cell.

implicit none
character*6, intent(in) :: side

type(auxi_crit_cell), pointer :: temp_n
type(geo_info) :: gn_cell
type(phy_info) :: pn_cell
type(state) :: sstt

nullify(temp_n)

! Find the stacked neighbor first.
select case(side)
 case('left  ')
  if(temp%l_stk%cv_nmb.gt.0) then
   call visit(temp%l_stk,temp_n)
   sstt=p_cell%l_state
  end if
 case('right ') 
  if(temp%r_stk%cv_nmb.gt.0) then
   call visit(temp%r_stk,temp_n)
   sstt=p_cell%r_state
  end if
 case default; call error_message
end select

! Update the hidden state in between and the ordinary cell-average in the 
! present and stacked discontinuity cells.
if(associated(temp_n)) then 
 gn_cell=temp_n%g_cell; pn_cell=temp_n%p_cell
 if(temp_n%l_stk.eq.address) then
  pn_cell%l_state=sstt
 else
  if(temp_n%r_stk.eq.address) then
   pn_cell%r_state=sstt
  else
   call error_message
  end if
 end if
 pn_cell%or_state=pn_cell%or_state+0.5d0*difference
 temp_n%p_cell=pn_cell
! Update the ordinary cell-average in the present critical cell.
 p_cell%or_state=p_cell%or_state+0.5d0*difference
end if

end subroutine update_stacked_neighbor


end subroutine return_data


subroutine update_xy(acv)
! Do the update in 'xy'-type critical cells, which is done by extrapolation of the updated 
! states in nearby 'xx'- and 'yy'-type critical cells.
implicit none
type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

integer :: i, j
type(geo_info) :: g_cell_p, g_cell_n
type(state), dimension(-3:3) :: uxy
type(state) :: x_state, y_state, updated_state
character*12 :: wave_type
character*14 :: fashion

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
 g_cell=temp%g_cell; p_cell=temp%p_cell; c_cell=temp%c_cell
 address=temp%address
 call clean_up_g_cell(g_cell_p); call clean_up_g_cell(g_cell_n)

 if(g_cell%g_type.eq.'xy ') then

! First decide in what fashion the updated-states should be taken.
  if(associated(temp%previous)) g_cell_p=temp%previous%g_cell
  if(associated(temp%next)) g_cell_n=temp%next%g_cell
  if(g_cell_p%g_type.eq.'xy '.and.g_cell_n%g_type.eq.'xy ') then
   fashion='no_direction'
  else
   if(g_cell_p%g_type.eq.'xy ') then
    if(g_cell%edge(2).eq.1.or.g_cell%edge(2).eq.3) then
     fashion='y_direction   '
    else
     fashion='x_direction   '
    end if
   else
    if(g_cell_n%g_type.eq.'xy ') then
     if(g_cell%edge(1).eq.1.or.g_cell%edge(1).eq.3) then
      fashion='y_direction   '
     else
      fashion='x_direction   '
     end if
    else
	 fashion='both_direction'
    end if
   endif
  end if
  if(g_cell_p%g_type.eq.'non'.and.g_cell_n%g_type.eq.'non') fashion='no_direction  '
  if(g_cell_p%g_type.eq.'xy '.and.g_cell_n%g_type.eq.'non') fashion='no_direction  '
  if(g_cell_p%g_type.eq.'non'.and.g_cell_n%g_type.eq.'xy ') fashion='no_direction  '
  
! Then evaluate the updated_states on the two sides.
  call evaluate_updated_state('left  ')
  c_cell%l_state(0)=updated_state
  call evaluate_updated_state('right ')
  c_cell%r_state(0)=updated_state

 end if

 temp%c_cell=c_cell

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do


contains


subroutine evaluate_updated_state(side)

implicit none
character*6, intent(in) :: side

character*10 :: x_direction, y_direction
type(state) :: l_state, r_state, c1, c2, l_phy, r_phy
real(8), dimension(2) :: normal
real(8) :: velocity, xxx, x_posi, y_posi, x_vel, y_vel, pressure
real(8) :: area_l, area_r
integer :: single
   
i=g_cell%x_idx; j=g_cell%y_idx
call find_single(g_cell,single)
   
select case(p_cell%wv_nb)
   
 case(1,3)
   
  select case(fashion)
   
   case('x_direction   ','y_direction   ','both_direction')
    select case(single)
     case(1)
      if(g_cell%point(single).eq.side) then
       x_direction='positive '; y_direction='positive '
      else
       x_direction='negative '; y_direction='negative '
      end if
     case(2)
      if(g_cell%point(single).eq.side) then
       x_direction='negative '; y_direction='positive '
      else
       x_direction='positive '; y_direction='negative '
      end if
     case(3)
      if(g_cell%point(single).eq.side) then
       x_direction='negative '; y_direction='negative '
      else
       x_direction='positive '; y_direction='positive '
      end if
     case(4)
      if(g_cell%point(single).eq.side) then
       x_direction='positive '; y_direction='negative '
      else
       x_direction='negative '; y_direction='positive '
      end if
     case default; call error_message
    end select    
  
    select case(fashion)
     case('both_direction')
      call smth_dpr(i,j,temp,uxy,'xx ',side, 'critical',.false.,3,'full    ')
      call inter_extrapolation(x_direction)       
      x_state=uxy(0)
      call smth_dpr(i,j,temp,uxy,'yy ',side, 'critical',.false.,3,'full    ')
      call inter_extrapolation(y_direction)       
      y_state=uxy(0)
      updated_state=0.5d0*(x_state+y_state)
     case('x_direction   ')
      call smth_dpr(i,j,temp,uxy,'xx ',side, 'critical',.false.,3,'full    ')
      call inter_extrapolation(x_direction)       
      x_state=uxy(0)
      updated_state=x_state
     case('y_direction   ')
      call smth_dpr(i,j,temp,uxy,'yy ',side, 'critical',.false.,3,'full    ')
      call inter_extrapolation(y_direction)       
      y_state=uxy(0)
      updated_state=y_state
    end select 
   case('no_direction  ')
    l_state=p_cell%l_state; r_state=p_cell%r_state
    normal=0.5d0*(g_cell%normal(1,:)+g_cell%normal(2,:))
     
    xxx=dsqrt(2.0d0)/2.0d0
    select case(single)
     case(1)
      if(g_cell%dis_posi(1).lt.-0.5d0.and.g_cell%dis_posi(2).lt.-0.5d0) then
       if(g_cell%edge(1).eq.1) then
        normal=(/xxx,xxx/)
       else
        normal=(/-xxx,-xxx/)
       end if
      end if
     case(2)
      if(g_cell%dis_posi(1).eq.1) then
       if(g_cell%dis_posi(1).gt.0.5d0.and.g_cell%dis_posi(2).lt.-0.5d0) then
        normal=(/xxx,-xxx/)
       end if
      else
       if(g_cell%dis_posi(1).lt.-0.5d0.and.g_cell%dis_posi(2).gt.0.5d0) then
        normal=(/-xxx,xxx/)
       end if
      end if
     case(3)
      if(g_cell%dis_posi(1).gt.0.5d0.and.g_cell%dis_posi(2).gt.0.5d0) then
       if(g_cell%edge(1).eq.2) then
        normal=(/xxx,xxx/)
       else
        normal=(/-xxx,-xxx/)
       end if
      end if
     case(4)
      if(g_cell%dis_posi(1).eq.1) then
       if(g_cell%dis_posi(1).gt.0.5d0.and.g_cell%dis_posi(2).lt.-0.5d0) then
        normal=(/xxx,-xxx/)
       end if
      else
       if(g_cell%dis_posi(1).lt.-0.5d0.and.g_cell%dis_posi(2).gt.0.5d0) then
        normal=(/-xxx,xxx/)
       end if
      end if
    end select
    
    x_posi=dfloat(g_cell%x_idx)*h; y_posi=dfloat(g_cell%y_idx)*h
    call riemann(l_state,r_state,normal,p_cell%wv_nb,x_posi,y_posi,c1,c2, &
                 velocity,wave_type)
    select case(side)
     case('left  '); updated_state=c1
     case('right '); updated_state=c2
     case default; call error_message
    end select
   
  end select  
  
 case(2)
    
  call find_physical_state(p_cell%l_state,(/1.0d0,0.0d0/),l_phy)
  call find_physical_state(p_cell%r_state,(/1.0d0,0.0d0/),r_phy)
  area_l=side_area(g_cell,'left  ')
  area_r=side_area(g_cell,'right ')
  x_vel=area_l*l_phy%value(2)+area_r*r_phy%value(2)
  y_vel=area_l*l_phy%value(3)+area_r*r_phy%value(3)
  pressure=area_l*l_phy%value(4)+area_r*r_phy%value(4)
   
  select case(side)
   case('left  ')
    l_phy%value(2)=x_vel
    l_phy%value(3)=y_vel
    l_phy%value(4)=pressure
    call find_conservative_state(l_phy,(/1.0d0,0.0d0/),updated_state)
   case('right ')		  
    r_phy%value(2)=x_vel
    r_phy%value(3)=y_vel
    r_phy%value(4)=pressure
    call find_conservative_state(r_phy,(/1.0d0,0.0d0/),updated_state)
   case default; call error_message
  end select
       
 case default; call error_message
    
end select
  
end subroutine evaluate_updated_state


subroutine inter_extrapolation(direction)

implicit none
character*10, intent(in) :: direction
logical :: ok
integer :: order

select case(direction)
 case('positive ')
  if(uxy(1).gt.0.9d0*error_data) then
   if(uxy(-2).gt.0.9d0*error_data) then
    uxy(0)=(-1.0d0/3.0d0)*uxy(-2)+uxy(-1)+1.0d0/3.0d0*uxy(1)
    order=2
   else
    if(uxy(-1).gt.0.9d0*error_data) then
     uxy(0)=0.5d0*(uxy(-1)+uxy(1))
     order=1
    else
     uxy(0)=uxy(1)
     order=0
    end if
   end if
   ok=whether_state_OK(uxy(0))
   if(.not.ok) then
    select case(order)
     case(2); uxy(0)=0.5d0*(uxy(-1)+uxy(1))
     case default; call error_message
    end select
    ok=whether_state_OK(uxy(0))
    if(.not.ok) call error_message
   end if
  else
   if(uxy(-3).gt.0.9d0*error_data) then
    uxy(0)=3.0d0*uxy(-1)-3.0d0*uxy(-2)+uxy(-3)
    order=2
   else
    if(uxy(-2).gt.0.9d0*error_data) then
     uxy(0)=2.0d0*uxy(-1)-uxy(-2)
     order=1
    else
     if(uxy(-1).gt.0.9d0*error_data) then
      uxy(0)=uxy(-1)
      order=0
     else
!      call error_message
     end if
    end if
   end if
   ok=whether_state_OK(uxy(0))
   if(.not.ok) then
    select case(order)
     case(2); uxy(0)=2.0d0*uxy(-1)-uxy(-2)
     case(1); uxy(0)=uxy(-1)
     case default; call error_message
    end select
    order=order-1
    ok=whether_state_OK(uxy(0))
    if(.not.ok) then
     select case(order)
      case(1); uxy(0)=uxy(-1)
      case default; call error_message
     end select
     order=order-1
     ok=whether_state_OK(uxy(0))
     if(.not.ok) call error_message
    end if     	  
   end if
  end if
  
 case('negative ')
  if(uxy(-1).gt.0.9d0*error_data) then
   if(uxy(2).gt.0.9d0*error_data) then
    uxy(0)=(-1.0d0/3.0d0)*uxy(2)+uxy(1)+1.0d0/3.0d0*uxy(-1)
    order=2
   else
    if(uxy(1).gt.0.9d0*error_data) then
     uxy(0)=0.5d0*(uxy(-1)+uxy(1))
     order=1
    else
     uxy(0)=uxy(-1)
     order=0
    end if
   end if
   ok=whether_state_OK(uxy(0))
   if(.not.ok) then
    select case(order)
     case(2); uxy(0)=0.5d0*(uxy(-1)+uxy(1))
     case default; call error_message
    end select
    ok=whether_state_OK(uxy(0))
    if(.not.ok) call error_message
   end if
  else
   if(uxy(3).gt.0.9d0*error_data) then
    uxy(0)=3.0d0*uxy(1)-3.0d0*uxy(2)+uxy(3)
    order=2
   else
    if(uxy(2).gt.0.9d0*error_data) then
     uxy(0)=2.0d0*uxy(1)-uxy(2)
     order=1
    else
     if(uxy(1).gt.0.9d0*error_data) then
      uxy(0)=uxy(1)
      order=0
     else
!      call error_message
     end if
    end if
   end if
   ok=whether_state_OK(uxy(0))
   if(.not.ok) then
    select case(order)
     case(2); uxy(0)=2.0d0*uxy(1)-uxy(2)
     case(1); uxy(0)=uxy(1)
     case default; call error_message
    end select
    order=order-1
    ok=whether_state_OK(uxy(0))
    if(.not.ok) then
     select case(order)
      case(1); uxy(0)=uxy(1)
      case default; call error_message
     end select
     order=order-1
     ok=whether_state_OK(uxy(0))
     if(.not.ok) call error_message
    end if     	  
   end if
  end if
  
 case default; call error_message

end select

end subroutine inter_extrapolation


end subroutine update_xy


end module update_left_n_right_states