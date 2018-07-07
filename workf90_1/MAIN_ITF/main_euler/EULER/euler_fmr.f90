! This module is for the physical state of the Euler equations. 

module physical_state

use euler_functions

use data_extrapolation
! 'data_ext.f90'

implicit none

type curve_map
 character*8 :: facing
 character*6 :: status
 integer :: number
end type

integer, parameter :: state_number=4, wave_number=3
! Parameters 'state_number' and 'wave_number' indicate the 
! possible numbers of physical states and waves. They may not
! be the same; for example, for the current Euler system in 2D, the 
! 'state_number' is 4 while 'wave_number' is 3. 

type state 
 real(8), dimension(4) :: value
end type state

interface assignment(=)
 module procedure assign_scalar, assign_euler
end interface

interface operator(.gt.)
 module procedure if_greater
end interface

interface operator(.lt.)
 module procedure if_less
end interface

interface operator(+)
 module procedure add
end interface

interface operator(-)
 module procedure subtract
end interface

interface operator(*)
 module procedure multiply, inner_product
end interface

interface operator(/)
 module procedure divide_r
end interface

interface weighted_division
 module procedure weighted_division_euler
end interface

interface data_fill
 module procedure data_fill_p
end interface

interface dabs
 module procedure absolute_value
end interface

public  f, g, state_number, wave_number, data_fill, riemann, &
        side_cut, find_physical_state
private add, subtract, multiply, divide_r, flux_in_x, flux_in_y, &
		alarm_in_density, triem, uf, vf, pf, cf, assign_scalar, &
		assign_euler, weighted_division_euler, inner_product 


contains


function f(u) result(c)
implicit none
type(state), intent(in) :: u
type(state) :: c
real(8), dimension(4) :: w
integer :: i
w=u%value
do i=1, 4
 c%value(i)=flux_in_x(i,w)
end do
end function f

function g(u) result(c)
implicit none
type(state), intent(in) :: u
type(state) :: c
real(8), dimension(4) :: w
integer :: i
w=u%value
do i=1, 4
 c%value(i)=flux_in_y(i,w)
end do
end function g


subroutine assign_scalar(a,value)

implicit none
type(state), intent(out) :: a
real(8), intent(in) :: value
integer :: i

do i=1, 4
 a%value(i)=value
end do
 
end subroutine assign_scalar


subroutine assign_euler(a,w)

implicit none
type(state), intent(out) :: a
type(state), intent(in) :: w

a%value=w%value
 
end subroutine assign_euler


function if_greater(a,value) result(c)

implicit none
type(state), intent(in) :: a
real(8), intent(in) :: value
logical :: c

integer :: i

c=.false.
do i=1,4 
 if(a%value(i).gt.value) then
  c=.true.; exit
 end if
end do

end function if_greater


function if_less(a,value) result(c)

implicit none
type(state), intent(in) :: a
real(8), intent(in) :: value
logical :: c

integer :: i

c=.false.
do i=1,4
 if(a%value(i).lt.value) then
  c=.true.
 end if
end do

end function if_less


function add(a,b) result(c)

implicit none
type(state), intent(in) :: a, b
type(state) :: c

c%value=a%value+b%value

end function add


function subtract(a,b) result(c)

implicit none
type(state), intent(in) :: a, b
type(state) :: c

c%value=a%value-b%value

end function subtract


function multiply(a,b) result(c)

implicit none
real(8), intent(in) :: a
type(state), intent(in) :: b
type(state) :: c

c%value=a*b%value

end function multiply


function inner_product(a,b) result(c)

implicit none
type(state), intent(in) :: a, b
real(8) :: c

integer :: i

c=0.0d0
do i=1, 4
 c=c+a%value(i)*b%value(i)
end do

end function inner_product


function divide_r(b,a) result(c)

implicit none
real(8), intent(in) :: a
type(state), intent(in) :: b
type(state) :: c

c%value=b%value/a

end function divide_r


function weighted_division_euler(yy,main_jump,left_difference,right_difference,normal, &
                                 wave_number) result(c)

implicit none
type(state), intent(in) :: yy, main_jump, left_difference, right_difference
real(8), dimension(2), intent(in) :: normal
integer, intent(in) :: wave_number
real(8) :: c

real(8), dimension(4) :: left_reliability, right_reliability, position
integer :: i
real(8), dimension(4) :: main_jump_normal_tangential
real(8), dimension(4) :: left_normal_tangential_difference
real(8), dimension(4) :: right_normal_tangential_difference
real(8), dimension(4) :: ordinary_normal_tangential
real(8) :: final_position, reliability, sum, aa, delta

! First, compute the vectors with normal and tangential components of momentums.
main_jump_normal_tangential(1)=main_jump%value(1)
main_jump_normal_tangential(2)=main_jump%value(2)*normal(1)+main_jump%value(3)*normal(2)
main_jump_normal_tangential(3)=-main_jump%value(2)*normal(2)+main_jump%value(3)*normal(1)
main_jump_normal_tangential(4)=main_jump%value(4)

left_normal_tangential_difference(1)=left_difference%value(1)
left_normal_tangential_difference(2)=left_difference%value(2)*normal(1)+ &
                                     left_difference%value(3)*normal(2)
left_normal_tangential_difference(3)=-left_difference%value(2)*normal(2)+ &
                                     left_difference%value(3)*normal(1)
left_normal_tangential_difference(4)=left_difference%value(4)

right_normal_tangential_difference(1)=right_difference%value(1)
right_normal_tangential_difference(2)=right_difference%value(2)*normal(1)+ &
                                     right_difference%value(3)*normal(2)
right_normal_tangential_difference(3)=-right_difference%value(2)*normal(2)+ &
                                     right_difference%value(3)*normal(1)
right_normal_tangential_difference(4)=right_difference%value(4)

ordinary_normal_tangential(1)=yy%value(1)
ordinary_normal_tangential(2)=yy%value(2)*normal(1)+yy%value(3)*normal(2)
ordinary_normal_tangential(3)=-yy%value(2)*normal(2)+yy%value(3)*normal(1)
ordinary_normal_tangential(4)=yy%value(4)

! Second, compute the left and right reliabilities.
do i=1, 4
 if(dabs(main_jump_normal_tangential(i)).lt.1.0d-4) then
  left_reliability(i)=1.0d0; right_reliability(i)=1.0d0
 else
  aa=left_normal_tangential_difference(i)/main_jump_normal_tangential(i)
  left_reliability(i)=2.0d0-dmin1(1.0d0,dabs(aa))
  aa=right_normal_tangential_difference(i)/main_jump_normal_tangential(i)
  right_reliability(i)=2.0d0-dmin1(1.0d0,dabs(aa))
 end if 
end do    

! Finally, compute the position
!sum=0.0d0; final_position=0.0d0
!do i=1,4
! reliability=dmin1(left_reliability(i),right_reliability(i))
! final_position=final_position+reliability*ordinary_normal_tangential(i)*main_jump_normal_tangential(i)
! sum=sum+reliability*main_jump_normal_tangential(i)*main_jump_normal_tangential(i)
!end do
!final_position=final_position/sum

select case(wave_number)
 case(1,3)
  sum=0.0d0; final_position=0.0d0
  do i=1,4
   reliability=dmin1(left_reliability(i),right_reliability(i))
   final_position=final_position+reliability*ordinary_normal_tangential(i)*main_jump_normal_tangential(i)
   sum=sum+reliability*main_jump_normal_tangential(i)*main_jump_normal_tangential(i)
  end do
  final_position=final_position/sum

 case(2)
  final_position=ordinary_normal_tangential(1)/main_jump_normal_tangential(1)

 case default; call error_message

end select   
c=final_position

end function weighted_division_euler


function absolute_value(a) result(c)

implicit none
type(state), intent(in) :: a
real(8) :: c

real(8), dimension(4) :: aa
integer :: i

do i=1,4
 aa(i)=dabs(a%value(i))
end do
c=dmax1(aa(1),aa(2),aa(3),aa(4))

end function absolute_value


function courant_number(a) result(c)

implicit none
type(state), intent(in) :: a
real(8) :: c

real(8), dimension(4) :: value
real(8) :: rho, u, v, cc

value=a%value
rho=value(1); u=uf(value); v=vf(value); cc=cf(value)
c=dmax1(dabs(u)+dabs(cc),dabs(v)+dabs(cc))

end function courant_number


subroutine riemann(wl,wr,normal,num,c1,c2,velocity,wave_type)
! Solution to Riemann problem on the direction of 'normal'.

implicit none
type(state), intent(in) :: wl, wr
real(8), dimension(2), intent(in) :: normal
integer, intent(in) :: num
type(state), intent(out) :: c1, c2
real(8), intent(out) :: velocity
character*12, intent(out) :: wave_type
! 'wl' and 'wr' are the states on the left and right side, 
! respectively, the 'normal' points from the right to the left.
! The output are the 'num'th and 'num+1'th states 'c1' and 'c2'
! in the Riemann separation, the 'wave_type', if it's a shock, a
! contact discontinuity, or a rarefaction wave, and the front 
! normal 'velocity'.

type(state) :: wml, wmr, mv
real(8) :: ss, p_front, p_rare
real(8), dimension(2) :: nnl
type(state) :: fl, fr, cc
integer :: i

real(8) :: rhoo, uu, vv, pp

do i=1,4
 c1%value(i)=-1000.0d0; c2%value(i)=-1000.0d0
 wml%value(i)=-1000.0d0; wmr%value(i)=-1000.0d0
end do
call triem(wl%value,wr%value,wml%value,wmr%value, &
           normal(1),normal(2))

 uu=uf(wml%value)
 vv=vf(wml%value)
 pp=pf(wml%value)

rhoo=wml%value(1)
uu=uf(wmr%value)
vv=vf(wmr%value)
pp=pf(wmr%value)
call tran3(rhoo,uu,vv,pp,cc%value)

! WL~WML LEFT-SHOCK, WML~CC SLIP, CC~WMR CONTACT DISCONTINUITY, WMR~WR RIGHT-SHOCK

select case(num)
 case(1); c1%value=wl%value; c2%value=wml%value
 case(2); c1%value=cc%value; c2%value=wmr%value
 case(3); c1%value=wmr%value; c2%value=wr%value
 case default; print*, 'Something must be wrong!'; pause
end select
 

select case(num)
 case(1)
  p_front=pf(c1%value); p_rare=pf(c2%value)
  if(p_front.gt.p_rare) then
   wave_type='rarefaction '
  else
   wave_type='shock       '
  end if
 case(2)
  wave_type='contact     '
 case(3)
  p_front=pf(c2%value); p_rare=pf(c1%value)
  if(p_front.gt.p_rare) then
   wave_type='rarefaction '
  else
   wave_type='shock       '
  end if
end select
  
nnl=-normal
ss=0.0d0
do i=1,4
 ss=ss+dabs(wl%value(i)-wr%value(i))
end do
select case(num)
 case(1,3)
  if(wave_type.eq.'shock       ') then
   if(ss.gt.0.1d-6) then
    fl=normal(1)*f(c1)+normal(2)*g(c1)
    fr=normal(1)*f(c2)+normal(2)*g(c2)
    velocity=(fl%value(1)-fr%value(1))/(c1%value(1)-c2%value(1))
   else
    mv=0.5d0*(c1+c2)
    velocity=normal(1)*uf(mv%value)+normal(2)*vf(mv%value)
    if(num.eq.1) then
     velocity=velocity-cf(mv%value) 
    else
     velocity=velocity+cf(mv%value)
    end if
   end if
  else
   velocity=-1000.0d0
  end if
 case(2)
  mv=0.5d0*(c1+c2)
  velocity=normal(1)*uf(mv%value)+normal(2)*vf(mv%value)
      	 
end select
 
end subroutine riemann


subroutine data_fill_p(uxy,highest_order)

implicit none

type(state), dimension(-3:3), intent(inout) :: uxy
integer, intent(in) :: highest_order

type(state), dimension(-3:3) :: uxy_temp
real(8), dimension(-3:3) :: uxyz
integer :: i, j, ho
logical :: if_ok 
real(8) :: psu

if_ok=.false.; ho=highest_order

! Extrapolate the physical states in a component-by-component way.
do while(.not.if_ok)
 do i=1, 4
  uxyz=error_data
  do j=-3,3 
   uxyz(j)=uxy(j)%value(i)
  end do
  call data_fill(uxyz,ho)
  do j=-3,3
   uxy_temp(j)%value(i)=uxyz(j)
  end do
 end do

! Check whether the extrapolated states are physical. If one of 
! them is not physical, the extrapolation must be redone with an
! order lower.
 if_ok=.true.
 do j=-ho,ho
  psu=pf(uxy_temp(j)%value)
  if(psu.le.1.0d-12.or.uxy_temp(j)%value(1).le.1.0d-12) then
   if_ok=.false.; exit
  end if
 end do
 if(.not.if_ok) ho=ho-1
 if(ho.le.0) then
  print*, 'Something must be wrong!'; pause
 end if

end do

do i=-3,3
 uxy(i)=uxy_temp(i)
end do
   
end subroutine data_fill_p


subroutine data_fill_flux(uxy,highest_order)

implicit none

type(state), dimension(-3:3), intent(inout) :: uxy
integer, intent(in) :: highest_order

real(8), dimension(-3:3) :: uxyz
integer :: i, j

do i=1, 4
 do j=-3,3 
  uxyz(j)=uxy(j)%value(i)
 end do
 call data_fill(uxyz,highest_order)
 do j=-3,3
  uxy(j)%value(i)=uxyz(j)
 end do
end do

end subroutine data_fill_flux


function wave_revers(wave) result(c)

implicit none
integer, intent(in) :: wave
integer :: c

c=4-wave

end function wave_revers


subroutine side_cut(state_1,state_2,state_or,difference)

implicit none

type(state), intent(in) :: state_1,state_2,state_or
type(state), intent(out) :: difference
integer :: i
logical, dimension(4) :: whether_cut

difference=0.0d0; whether_cut=.false.

do i=1,4
 if(state_1%value(i).lt.state_2%value(i)) then
  if(state_or%value(i).gt.state_2%value(i)) then
   whether_cut(i)=.true.
  end if
 else
  if(state_1%value(i).gt.state_2%value(i)) then
   if(state_or%value(i).lt.state_2%value(i)) then
    whether_cut(i)=.true.
   end if
  end if
 end if
end do

do i=1,4
 if(whether_cut(i)) difference%value(i)=state_or%value(i)-state_2%value(i)
end do

end subroutine side_cut


function whether_state_OK(u) result(c)
! This function checks whether the given state has positive density, energy and pressure.

implicit none
type(state), intent(in) :: u
logical :: c

real(8), dimension(4) :: value
real(8) :: p
integer :: i

do i=1,4
 value(i)=u%value(i)
end do
c=.true.
if(value(1).le.0.0d0) then
 c=.false.; return
end if
if(value(4).le.0.0d0) then
 c=.false.; return
end if
p=pf(value)
if(p.le.0.0d0) then
 c=.false.; return
end if

end function whether_state_OK


subroutine find_physical_state(u,normal,uu)
! Given a conservation state, this subroutine gives out the density, 
! x- and y-velocities and pressure.

implicit none

type(state), intent(in) :: u
real(8), dimension(2), intent(in) :: normal
type(state), intent(out) :: uu

real(8), dimension(4) :: v, vv

v=u%value
vv(1)=v(1)
vv(2)=normal(1)*uf(v)+normal(2)*vf(v)
vv(3)=-normal(2)*uf(v)+normal(1)*vf(v)
vv(4)=pf(v)
uu%value=vv

end subroutine find_physical_state


subroutine find_conservative_state(uu,normal,u)
! Given a physical state in density, x- and y-velocities and energy, this subroutine gives
! out the corresponding conservative state.

implicit none

type(state), intent(in) :: uu
real(8), dimension(2), intent(in) :: normal
type(state), intent(out) :: u

real(8), dimension(4) :: vv
real(8) :: rhoi, ui, vi, pi

vv=uu%value
rhoi=vv(1)
ui=normal(1)*vv(2)-normal(2)*vv(3)
vi=normal(2)*vv(2)+normal(1)*vv(3)
pi=vv(4)
call tran3(rhoi,ui,vi,pi,vv)
u%value=vv


end subroutine find_conservative_state


end module physical_state