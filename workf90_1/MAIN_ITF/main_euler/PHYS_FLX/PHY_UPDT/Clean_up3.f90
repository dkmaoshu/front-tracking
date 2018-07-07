module clean_up
! This module is for the clean_up procedure of the discontinuity curves.

use discontinuity_curves
! 'discontinuity_curve.f90'

implicit none

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(state) :: l_store, r_store

public  clean_up_discontinuities
private g_cell, p_cell, l_store, r_store, cleaning_data, return_data, &
        uniform_states_between


contains


subroutine clean_up_discontinuities
! Implement the clean_up.

implicit none

integer :: i

do i=1, cvn
 if(cvv(i)%status.eq.'awake ') then
  call uniform_states_between(i)
 end if
end do

! call check_list_c(1,'down  ')

do i=1, cvn
 if(cvv(i)%status.eq.'awake ') then
  call cleaning_data(i)
 end if
end do

! call check_list_c(1,'down  ')

do i=1, cvn
 if(cvv(i)%status.eq.'awake ') then
  call return_data(i)
 end if
end do

! call check_list_c(1,'down  ')


end subroutine clean_up_discontinuities


subroutine uniform_states_between(cv_nb)

implicit none

integer, intent(in) :: cv_nb

type(adss_info) :: address
integer :: head_mark
logical :: head_switch
type(critical_cell), pointer :: temp_n
type(state) :: state_in_between
type(geo_info) :: g_cell

nullify(temp)
if(associated(cvv(cv_nb)%begin%next)) then
 temp=>cvv(cv_nb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch) 
 
 g_cell=temp%g_cell; address=temp%address
 if(temp%l_stk.ne.adss_info(0,0)) then
  call visit(temp%l_stk,temp_n)
  if(temp_n%l_stk.eq.address) then
   state_in_between=(temp%p_cell%l_state+temp_n%p_cell%l_state)/2.0d0
   temp_n%p_cell%l_state=state_in_between
  else
   if(temp_n%r_stk.eq.address) then
    state_in_between=(temp%p_cell%l_state+temp_n%p_cell%r_state)/2.0d0
    temp_n%p_cell%r_state=state_in_between
   else
    call error_message
   end if
  end if
  temp%p_cell%l_state=state_in_between
 end if

 if(temp%r_stk.ne.adss_info(0,0)) then
  call visit(temp%r_stk,temp_n)
  if(temp_n%l_stk.eq.address) then
   state_in_between=(temp%p_cell%r_state+temp_n%p_cell%l_state)/2.0d0
   temp_n%p_cell%l_state=state_in_between
  else
   if(temp_n%r_stk.eq.address) then
    state_in_between=(temp%p_cell%r_state+temp_n%p_cell%r_state)/2.0d0
    temp_n%p_cell%r_state=state_in_between
   else
    call error_message
   end if
  end if  
  temp%p_cell%r_state=state_in_between
 end if

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine uniform_states_between


subroutine cleaning_data(cv_nmb)

implicit none
integer, intent(in) :: cv_nmb

type(adss_info) :: address
integer :: head_mark
logical :: head_switch

type(geo_info) :: g_cell_n
type(state) :: l_state, r_state
type(state) :: left_phy_state, right_phy_state, ordin_state, ordin_phy_state
type(state) :: difference, difference_normal, left_difference, right_difference
type(state) :: left_difference_normal, right_difference_normal
real(8) :: l_x_velo, r_x_velo, or_x_velo, l_y_velo, r_y_velo, or_y_velo, l_pressure, &
           r_pressure, or_pressure, or_density, ttt
real(8), dimension(2) :: normal 
real(8) :: l_area, r_area, velocity, alpha, beta
character*12 :: wave_type

type(state), dimension(4) :: vec, lleft_eigenvectors, lright_eigenvectors,  &
                             rleft_eigenvectors, rright_eigenvectors
type(state) :: ln, rn, c1, c2, nc1, nc2, dd
integer :: i, j

real(8) :: x_posi, y_posi

type(critical_cell), pointer :: temp_m
type(state) :: ss_phy_state
real(8) :: inteng_l, inteng_r, inteng_m, inteng_or, inteng_s, m_area, m_pressure

nullify(temp)
if(associated(cvv(cv_nmb)%begin%next)) then
 temp=>cvv(cv_nmb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

difference=error_data; difference_normal=error_data
left_phy_state=error_data; right_phy_state=error_data; ordin_phy_state=error_data
ordin_state=error_data
left_difference=error_data; right_difference=error_data
left_difference_normal=error_data; right_difference_normal=error_data
l_x_velo=error_data; r_x_velo=error_data; or_x_velo=error_data
l_y_velo=error_data; r_y_velo=error_data; or_y_velo=error_data
l_pressure=error_data; r_pressure=error_data; or_pressure=error_data
normal=error_data
ln=error_data; rn=error_data; c1=error_data; c2=error_data
nc1=error_data; nc2=error_data
dd=error_data

do i=1, 4
 vec(i)=error_data
 lleft_eigenvectors(i)=error_data
 lright_eigenvectors(i)=error_data
 rleft_eigenvectors(i)=error_data
 rright_eigenvectors(i)=error_data
end do

do while(associated(temp).and.head_switch) 
 g_cell=temp%g_cell; p_cell=temp%p_cell
 address=temp%address
 x_posi=dfloat(g_cell%x_idx)*h; y_posi=dfloat(g_cell%y_idx)*h

 left_difference=error_data
 left_difference_normal=error_data
 right_difference=error_data
 right_difference_normal=error_data
 difference=error_data
 difference_normal=error_data
  
 select case(p_cell%wv_nb)
! The cleaning procedures for shocks and contact discontinuities
! are different.
  
  case(1,3)
   
! First, compute the normal of the critical cell.    
   call compute_normal(normal)
    
! Then compute the conservation difference.
   l_state=p_cell%l_state; l_store=l_state
   r_state=p_cell%r_state; r_store=r_state
   l_area=side_area(g_cell,'left  ')
   r_area=side_area(g_cell,'right ')
   difference=l_area*l_state+r_area*r_state
   difference=p_cell%or_state-difference
   
! And then decompose the conservation difference according to the characteristics.
   call transform_1(l_state,normal,ln)
   call transform_1(r_state,normal,rn)
   call transform_1(difference,normal,difference_normal)
    
   select case(p_cell%wv_nb)
    
    case(1)
     call roe_decomposition(rn,rn,rleft_eigenvectors,rright_eigenvectors)
     do i=2, 4
      vec(i)=rright_eigenvectors(i)
     end do
      
     call riemann(l_state,r_state,normal,1,x_posi,y_posi,c1,c2,velocity,wave_type)
     call transform_1(c1,normal,nc1)
     call transform_1(c2,normal,nc2)
     vec(1)=nc2-nc1
     alpha=inner_product_state(rleft_eigenvectors(1),difference_normal)
     alpha=alpha/inner_product_state(vec(1),rleft_eigenvectors(1))
     right_difference_normal=difference_normal-alpha*vec(1)
     dd=0.0d0
     do i=2,4
      beta=inner_product_state(right_difference_normal,rleft_eigenvectors(i))
      beta=beta/inner_product_state(vec(i),lleft_eigenvectors(i))
!      beta=dsign(1.0d0,beta)*dmin1(dabs(beta),dabs(beta)**(3.0d0/2.0d0))
      dd=dd+beta*vec(i)
     end do
     right_difference_normal=dd    
     call transform_2(right_difference_normal,normal,right_difference)
     r_store=r_store+right_difference
     
    case(3) 
     call roe_decomposition(ln,ln,lleft_eigenvectors,lright_eigenvectors)
     do i=1, 3
      vec(i)=lright_eigenvectors(i)
     end do
     call riemann(l_state,r_state,normal,3,x_posi,y_posi,c1,c2,velocity,wave_type)
     call transform_1(c1,normal,nc1)
     call transform_1(c2,normal,nc2)
     vec(4)=nc2-nc1
     alpha=inner_product_state(lleft_eigenvectors(4),difference_normal)
     alpha=alpha/inner_product_state(vec(4),lleft_eigenvectors(4))
     left_difference_normal=difference_normal-alpha*vec(4)
     dd=0.0d0
     do i=1,3
      beta=inner_product_state(left_difference_normal,lleft_eigenvectors(i))
      beta=beta/inner_product_state(vec(i),lleft_eigenvectors(i))
!      beta=dsign(1.0d0,beta)*dmin1(dabs(beta),dabs(beta)**(3.0d0/2.0d0))
      dd=dd+beta*vec(i)
     end do
     left_difference_normal=dd    
     call transform_2(left_difference_normal,normal,left_difference)
     l_store=l_store+left_difference
     
    end select
                  
  case(2)
   
! Slip lines are ignored in the present version of code. 

   ordin_state=p_cell%or_state

! For stacked discontinuity cells the velocity and pressure is uniform. To
! compute them, we need to compute the ordinary cell-average over the 
! stacked grid cell.

   if(temp%l_stk.ne.adss_info(0,0)) call accumulate('left  ')
   if(temp%r_stk.ne.adss_info(0,0)) call accumulate('right ')

   call find_physical_state(p_cell%l_state,(/1.0d0,0.0d0/),left_phy_state)
!   l_x_velo=left_phy_state%value(2)
!   l_y_velo=left_phy_state%value(3)
!   l_pressure=left_phy_state%value(4)
   call find_physical_state(p_cell%r_state,(/1.0d0,0.0d0/),right_phy_state)
!   r_x_velo=right_phy_state%value(2)
!   r_y_velo=right_phy_state%value(3)
!   r_pressure=right_phy_state%value(4)
!   call find_physical_state(ordin_state,(/1.0d0,0.0d0/),ordin_phy_state)
!   or_x_velo=ordin_phy_state%value(2)
!   or_y_velo=ordin_phy_state%value(3)
!   or_pressure=ordin_phy_state%value(4) 

   or_x_velo=x_velocity(ordin_state)
   or_y_velo=y_velocity(ordin_state)
   inteng_or=internal_energy(ordin_state)
   
!  important changes for two different gases 
      
   l_area=side_area(g_cell,'left  ')
   r_area=1.0d0-l_area
   inteng_m=0.0d0; m_area=0.0d0; m_pressure=0.0d0 

   inteng_l=internal_energy(p_cell%l_state)
   l_pressure=pressure(p_cell%l_state)

   inteng_r=internal_energy(p_cell%r_state)
   r_pressure=pressure(p_cell%r_state)

!   inteng_or=internal_energy(p_cell%or_state)
    
   if(temp%l_stk.ne.adss_info(0,0)) then
    inteng_m=inteng_l; m_pressure=l_pressure; m_area=l_area
    call visit(temp%l_stk,temp_m) 
    if(temp_m%l_stk.eq.temp%address) then
     l_area=side_area(temp_m%g_cell,'right ')
     inteng_l=internal_energy(temp_m%p_cell%r_state)
     l_pressure=pressure(temp_m%p_cell%r_state)
     m_area=m_area-l_area
    else
     if(temp_m%r_stk.eq.temp%address) then
      l_area=side_area(temp_m%g_cell,'left  ')
      inteng_l=internal_energy(temp_m%p_cell%l_state)
	  l_pressure=pressure(temp_m%p_cell%l_state)
      m_area=m_area-l_area 
     else
      call error_message
     end if
 !    inteng_or=inteng_or+internal_energy(temp_m%p_cell%or_state)-inteng_m
    end if
   end if

   if(temp%r_stk.ne.adss_info(0,0)) then
    inteng_m=inteng_r; m_pressure=r_pressure; m_area=r_area
    call visit(temp%r_stk,temp_m) 
    if(temp_m%l_stk.eq.temp%address) then
     r_area=side_area(temp_m%g_cell,'right ')
     inteng_r=internal_energy(temp_m%p_cell%r_state)
	 r_pressure=pressure(temp_m%p_cell%r_state)
     m_area=m_area-r_area
    else
     if(temp_m%r_stk.eq.temp%address) then
      r_area=side_area(temp_m%g_cell,'left  ')
      inteng_r=internal_energy(temp_m%p_cell%l_state)
	  r_pressure=pressure(temp_m%p_cell%l_state)
      m_area=m_area-r_area
     else
      call error_message
     end if
!     inteng_or=inteng_or+internal_energy(temp_m%p_cell%or_state)-inteng_m
    end if
   end if  

   or_pressure=l_area*l_pressure+m_area*m_pressure+r_area*r_pressure
   inteng_s=l_area*inteng_l+m_area*inteng_m+r_area*inteng_r
   or_pressure=or_pressure*inteng_or/inteng_s   
           
!   or_density=l_area*p_cell%l_state%value(1)+r_area*p_cell%r_state%value(1)
!   or_pressure=or_pressure*ordin_state%value(1)/or_density

   ttt=or_x_velo-l_x_velo
!   ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),dabs(ttt)**(3.0d0/2.0d0))
   l_x_velo=l_x_velo+ttt
   ttt=or_x_velo-r_x_velo
!   ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),dabs(ttt)**(3.0d0/2.0d0))
   r_x_velo=r_x_velo+ttt
  
   ttt=or_y_velo-l_y_velo
!   ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),dabs(ttt)**(3.0d0/2.0d0))
   l_y_velo=l_y_velo+ttt
   ttt=or_y_velo-r_y_velo
!   ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),dabs(ttt)**(3.0d0/2.0d0))
   r_y_velo=r_y_velo+ttt
   
   ttt=or_pressure-l_pressure
!   ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),dabs(ttt)**(3.0d0/2.0d0))
   l_pressure=l_pressure+ttt
   ttt=or_pressure-r_pressure
!   ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),dabs(ttt)**(3.0d0/2.0d0))
   r_pressure=r_pressure+ttt
  
   left_phy_state%value(2)=l_x_velo
   left_phy_state%value(3)=l_y_velo
   left_phy_state%value(4)=l_pressure
   right_phy_state%value(2)=r_x_velo
   right_phy_state%value(3)=r_y_velo
   right_phy_state%value(4)=r_pressure
  
   call find_conservative_state(left_phy_state,(/1.0d0,0.0d0/),l_store)
   call find_conservative_state(right_phy_state,(/1.0d0,0.0d0/),r_store)
  
 end select
 
 temp%l_store=l_store; temp%r_store=r_store

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do


contains


subroutine compute_normal(normal)

implicit none

real(8), dimension(2), intent(out) :: normal(2)

type(geo_info) :: g_cell_p, g_cell_n
real(8), dimension(2) :: pt1, pt2
real(8) :: rr
logical :: previous_not_end, next_not_end
integer :: single

previous_not_end=.false.; next_not_end=.false.

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)
rr=length_of_segment(pt1,pt2)
if(rr.lt.1.0d-2) then
 if(associated(temp%previous)) then
  g_cell_p=temp%previous%g_cell
  call pick_point(g_cell_p,1,pt1)
  previous_not_end=.true. 
 end if
 if(associated(temp%next)) then
  g_cell_n=temp%next%g_cell
  call pick_point(g_cell_n,2,pt2)
 end if
 if(.not.previous_not_end.and..not.next_not_end) then
  if(g_cell%g_type.ne.'xy ') call error_message
  call find_single(g_cell,single)
  select case(single)
   case(1)
    if(g_cell%edge(1).eq.1) then
     pt1=(/1.0d0,0.0d0/); pt2=(/0.0d0,1.0d0/) 
    else
     pt1=(/0.0d0,1.0d0/); pt2=(/1.0d0,0.0d0/)
    end if
   case(2)
    if(g_cell%edge(1).eq.1) then
     pt1=(/-1.0d0,0.0d0/); pt2=(/0.0d0,1.0d0/) 
    else
     pt1=(/0.0d0,1.0d0/); pt2=(/-1.0d0,0.0d0/)
    end if
   case(3)
    if(g_cell%edge(1).eq.3) then
     pt1=(/-1.0d0,0.0d0/); pt2=(/0.0d0,-1.0d0/) 
    else
     pt1=(/0.0d0,-1.0d0/); pt2=(/-1.0d0,0.0d0/)
    end if
   case(4)
    if(g_cell%edge(1).eq.3) then
     pt1=(/1.0d0,0.0d0/); pt2=(/0.0d0,-1.0d0/) 
    else
     pt1=(/0.0d0,-1.0d0/); pt2=(/1.0d0,0.0d0/)
    end if
  end select 	 
 end if
end if

normal=normal_of_line(pt1,pt2)

end subroutine compute_normal


subroutine transform_1(given_state,normal,transformed_state)
! This subroutine transforms the 'given_state' in 'xy'-coordinates into 
! the 'transformed_state' in 'normal-tangential'-coordinates.

implicit none
type(state), intent(in) :: given_state
real(8), dimension(2) :: normal
type(state), intent(out) :: transformed_state

real(8), dimension(2) :: mxy, mnt, tangent

tangent(1)=-normal(2); tangent(2)=normal(1)
transformed_state=given_state
mxy(1)=given_state%value(2); mxy(2)=given_state%value(3)
mnt(1)=inner_product(normal,mxy)
mnt(2)=inner_product(tangent,mxy)
transformed_state%value(2)=mnt(1)
transformed_state%value(3)=mnt(2)

end subroutine transform_1


subroutine transform_2(given_state,normal,transformed_state)
! This subroutine transforms the 'given_state' in 'xy'-coordinates into
! the 'transformed_state' in 'normal-tangential' coordinates.

implicit none
type(state), intent(in) :: given_state
real(8), dimension(2) :: normal
type(state), intent(out) :: transformed_state

real(8), dimension(2) :: mxy, mnt, tangent

tangent(1)=-normal(2); tangent(2)=normal(1)
transformed_state=given_state
mnt(1)=given_state%value(2)
mnt(2)=given_state%value(3)
mxy(1)=mnt(1)*normal(1)+mnt(2)*tangent(1)
mxy(2)=mnt(1)*normal(2)+mnt(2)*tangent(2)
transformed_state%value(2)=mxy(1)
transformed_state%value(3)=mxy(2)

end subroutine transform_2


function inner_product(vec1,vec2) result(c)

implicit none
real(8), dimension(:), intent(in) :: vec1, vec2
real(8) :: c

integer :: size1, size2, i

size1=size(vec1); size2=size(vec2)
if(size1.ne.size2) call error_message

c=0.0d0
do i=1,size1
 c=c+vec1(i)*vec2(i)
end do

end function inner_product


function inner_product_state(state_1,state_2) result(c)

implicit none

type(state), intent(in) :: state_1, state_2
real(8) :: c

real(8), dimension(4) :: w1, w2
w1=state_1%value; w2=state_2%value
c=inner_product(w1,w2)

end function inner_product_state


subroutine roe_decomposition(state_1,state_2,left_eigenvectors,right_eigenvectors)
! Given physical states 'state_1' and 'state_2', finds the eigenvectors of
! the corresponding Roe's Riemann approximation solver.

implicit none
type(state), intent(in) :: state_1, state_2
type(state), dimension(4), intent(out) :: left_eigenvectors, right_eigenvectors

real(8), dimension(4) :: w1, w2, ch, coef
real(8), dimension(4,4) :: el, er
integer :: i

w1=state_1%value; w2=state_2%value
call roe(w2,w1,el,er,ch,coef,state_1%gamma)
do i=1,4
 left_eigenvectors(i)%value=el(i,:)
 right_eigenvectors(i)%value=er(i,:)
end do

end subroutine roe_decomposition


subroutine accumulate(side)

implicit none

character*6, intent(in) :: side

type(critical_cell), pointer :: temp_n, temp_nn
character*6 :: side_local
logical :: go_on

side_local=side
go_on=.true.; temp_n=>temp; nullify(temp_nn)
do while(go_on) 
 select case(side_local)
  case('left  ') 
   call visit(temp_n%l_stk,temp_nn)
   if(temp_nn%p_cell%wv_nb.eq.2) then
    ordin_state=ordin_state+temp_nn%p_cell%or_state-temp_n%p_cell%l_state
   end if 
  case('right ') 
   call visit(temp_n%r_stk,temp_nn)
   if(temp_nn%p_cell%wv_nb.eq.2) then 
    ordin_state=ordin_state+temp_nn%p_cell%or_state-temp_n%p_cell%r_state
   end if
  case default; call error_message
 end select
 if(temp_nn%p_cell%wv_nb.ne.2) return
 if(temp_nn%l_stk.eq.temp_n%address) then
  if(temp_nn%r_stk.ne.adss_info(0,0)) then
   temp_n=>temp_nn; nullify(temp_nn); side_local='right '
  else
   return
  end if
 else
  if(temp_nn%r_stk.ne.temp_n%address) call error_message
  if(temp_nn%l_stk.ne.adss_info(0,0)) then
   temp_n=>temp_nn; nullify(temp_nn); side_local='left  '
  else
   return
  end if
 end if
end do

end subroutine accumulate


end subroutine cleaning_data


subroutine return_data(cv_nb)
! Return the data stored in 'c_cell' back to the left and right states in 'xx'- and 'yy'-type
! critical cells.
implicit none
integer, intent(in) :: cv_nb

type(adss_info) :: address !, address1
integer :: head_mark
logical :: head_switch
type(state) :: l_store, r_store, difference, or_state, phy_st

! real(8) :: xx, yy, rr
! real(8) :: u_radial_l, u_radial_r, u_tangen_l, u_tangen_r

nullify(temp)
if(associated(cvv(cv_nb)%begin%next)) then
 temp=>cvv(cv_nb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch) 
 g_cell=temp%g_cell; p_cell=temp%p_cell
 l_store=temp%l_store; r_store=temp%r_store
 address=temp%address
  
 difference=error_data; or_state=error_data
 select case(p_cell%wv_nb)
   
  case(1,3)
   difference=l_store-p_cell%l_state
   p_cell%l_state=l_store
   
  case(2)
!   call find_physical_state(l_store,(/1.0d0,0.0d0/),phy_st)
!   phy_st%value(1)=p_cell%or_state%value(1)
!   call find_conservative_state(phy_st,(/1.0d0,0.0d0/),or_state)
!   difference=(p_cell%or_state-p_cell%l_state)-(or_state-l_store)
   difference=l_store-p_cell%l_state
   p_cell%l_state=l_store
!   p_cell%or_state=or_state
   
  case default; call error_message
   
 end select
  
 call update_stacked_neighbor('left  ')
   
 difference=error_data; or_state=error_data
 select case(p_cell%wv_nb)
   
  case(1,3)
   difference=r_store-p_cell%r_state
   p_cell%r_state=r_store
   
  case(2)
!   call find_physical_state(r_store,(/1.0d0,0.0d0/),phy_st)
!   phy_st%value(1)=p_cell%or_state%value(1)
!   call find_conservative_state(phy_st,(/1.0d0,0.0d0/),or_state)
!   difference=(p_cell%or_state-p_cell%r_state)-(or_state-r_store)
   difference=r_store-p_cell%r_state
   p_cell%r_state=r_store
!   p_cell%or_state=or_state
  
  case default; call error_message
   
 end select

 call update_stacked_neighbor('right ')
  
 temp%p_cell=p_cell
  
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do


contains


subroutine update_stacked_neighbor(side)
! This subroutine updates the left or right states in neighboring critical cells which are
! stacked in the same grid cell. It also updates the ordinary cell-average in the present 
! critical cell.

implicit none
character*6, intent(in) :: side

type(critical_cell), pointer :: temp_n
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

! Update the hidden state in between and the ordinary cell-average in 
! the present critical cell.
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

 select case(p_cell%wv_nb)
  case(1,3)
   pn_cell%or_state=pn_cell%or_state+0.5d0*difference
   temp_n%p_cell=pn_cell
! Update the ordinary cell-average in the present critical cell.
   p_cell%or_state=p_cell%or_state+0.5d0*difference
  case(2)
   select case(pn_cell%wv_nb)
    case(1,3)
     pn_cell%or_state=pn_cell%or_state+difference
    case(2)
    case default; call error_message
   end select
  case default; call error_message
 end select
end if
  
end subroutine update_stacked_neighbor


end subroutine return_data


end module clean_up