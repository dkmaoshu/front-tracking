module evaluate_updated_xxyy

use produce_polynomial_4_xx_yy
! 'proc_pln.f90'

use euler_functions
! 'euler_functions.f90'

type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(comp_info) :: c_cell

real(8), dimension(-1:1) :: x, y
type(state), dimension(-1:1) :: left_nearest, right_nearest
type(state) :: l_state, r_state, l_updated_state, r_updated_state
real(8), dimension(-1:1,-1:1,2) :: integral_points
type(state), dimension(-1:1,-1:1) :: integral_data
real(8), dimension(0:2) :: interpolation_points
type(state), dimension(0:2) :: interpolation_data
real(8), dimension(-1:1,2) :: normal

public  update_xx_yy
private	x, y, left_nearest, right_nearest, l_state, r_state, &
        updated_state, integral_points, integral_data, &
		interpolation_points, interpolation_data,  &
		evaluate_updated_state, simpson_integral, &
		produce_interpolation_points, produce_interpolation_data, &
		produce_integral_points, integral_state, normal, &
		l_updated_state, r_updated_state


contains


subroutine update_xx_yy(acv)
! Do the update in 'xx'- and 'yy'-type critical cells, which is 
! done by extrapolation and numerical integral.

implicit none
type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark
logical :: head_switch

integer :: i
real(8), dimension(2) :: pt, nml
real(8) :: velocity, rrr 
character*12 :: wave_type

! type(geo_info) :: gggxxx

nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch) 
 g_cell=temp%g_cell; p_cell=temp%p_cell; c_cell=temp%c_cell
 address=temp%address

!  gggxxx=g_cell
 
 if(g_cell%g_type.eq.'xx '.or.g_cell%g_type.eq.'yy ') then
  
! Compute the nearest states on the curve segment.
  call pick_point(g_cell,1,pt)
  x(-1)=pt(1); y(-1)=pt(2)
  normal(-1,:)=g_cell%normal(1,:)
  call pick_point(g_cell,2,pt)
  x(1)=pt(1); y(1)=pt(2)
  normal(1,:)=g_cell%normal(2,:)
  x(0)=g_cell%mdis_posi(1); y(0)=g_cell%mdis_posi(2)
  normal(0,:)=0.5d0*(normal(-1,:)+normal(1,:))
  rrr=normal(0,1)*normal(0,1)+normal(0,2)*normal(0,2)
  rrr=dsqrt(rrr)
  normal(0,:)=normal(0,:)/rrr   
  do i=-1,1
   l_state=p_side(temp,x(i),y(i),'left  ')
   r_state=p_side(temp,x(i),y(i),'right ')
   nml=normal(i,:)
   call riemann(l_state,r_state,nml,p_cell%wv_nb,left_nearest(i),right_nearest(i),velocity, &
               wave_type)
  end do
 
! Evaluate the updated reconstruction solution in the critical cell and store the data in
! 'c_cell'.
 call evaluate_updated_state

 end if

 temp%c_cell=c_cell

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine update_xx_yy


subroutine evaluate_updated_state
! It evaluates the updated left and right states. It consists of 
! two steps: integration and clean-up.

implicit none

type(state) :: left_integral, right_integral, l_original_state, &
               r_original_state, lno, lnu, rno, rnu, difference, difference_normal
type(state) :: left_difference, right_difference, left_difference_normal, &
			   right_difference_normal 
real(8), dimension(2) :: nml
type(state), dimension(4) :: vec, lleft_eigenvectors, lright_eigenvectors,  &
                             rleft_eigenvectors, rright_eigenvectors

 type(state), dimension(4) :: lm, rm
!  type(state) :: xxx, yyy
 real(8) :: checkk

type(state) :: left_nearest_in_normal, right_nearest_in_normal,  &
               intermediate_nearest_in_normal, dd

real(8) :: alpha, beta, limiter, delta
real(8), dimension(3,3) :: coefficients
real(8), dimension(3) :: strength, rhs, sma
integer :: i, j

type(state) :: left_phy_state, right_phy_state, ordin_state, ordin_phy_state
real(8) :: l_tangential_velo, r_tangential_velo, or_tangential_velo
real(8) :: l_normal_velo, r_normal_velo, or_normal_velo
real(8) :: l_pressure, r_pressure, or_pressure, ttt, position


! type(geo_info) :: check
! type(state) :: xxx
  
! check=g_cell

left_difference=error_data
left_difference_normal=error_data
right_difference=error_data
right_difference_normal=error_data

! First, produce tentative updated states.
call produce_integral_points('both  ')
l_original_state=p_cell%l_state 
call integral_state('left  ','updated ',l_updated_state)
r_original_state=p_cell%r_state
call integral_state('right ','updated ',r_updated_state)

nml=normal(0,:)
call transform_1(l_original_state,nml,lno)
call transform_1(l_updated_state,nml,lnu)
call transform_1(r_original_state,nml,rno)
call transform_1(r_updated_state,nml,rnu)

call produce_integral_points('left  ')
call integral_state('left  ','updated ',left_integral)
call produce_integral_points('right ')
call integral_state('right ','updated ',right_integral)
difference=p_cell%or_state-left_integral-right_integral
call transform_1(difference,nml,difference_normal)

! Computer the limiter.
call produce_limiter

select case(p_cell%wv_nb)
 case(1)
  call roe_decomposition(rno,rnu,rleft_eigenvectors,rright_eigenvectors)
  do i=2, 4
   vec(i)=rright_eigenvectors(i)
  end do
  call transform_1(left_nearest(0),nml,left_nearest_in_normal)
  call transform_1(right_nearest(0),nml,right_nearest_in_normal)
  vec(1)=right_nearest_in_normal-left_nearest_in_normal
  alpha=inner_product_state(rleft_eigenvectors(1),difference_normal)
  alpha=alpha/inner_product_state(vec(1),rleft_eigenvectors(1))
  right_difference_normal=difference_normal-alpha*vec(1)
  dd=0.0d0
  do i=2,4
   beta=inner_product_state(right_difference_normal, &
                               lleft_eigenvectors(i))
   beta=beta/inner_product_state(vec(i),lleft_eigenvectors(i))
   beta=dsign(1.0d0,beta)*dmin1(dabs(beta),limiter)
   dd=dd+beta*vec(i)
  end do
  right_difference_normal=dd    
  call transform_2(right_difference_normal,nml,right_difference)
  r_updated_state=r_updated_state+right_difference

 case(3)
  call roe_decomposition(lno,lnu,lleft_eigenvectors,lright_eigenvectors)
  do i=1, 3
   vec(i)=lright_eigenvectors(i)
  end do
  call transform_1(left_nearest(0),nml,left_nearest_in_normal)
  call transform_1(right_nearest(0),nml,right_nearest_in_normal)
  vec(4)=right_nearest_in_normal-left_nearest_in_normal
  alpha=inner_product_state(lleft_eigenvectors(4),difference_normal)
  alpha=alpha/inner_product_state(vec(4),lleft_eigenvectors(4))
  left_difference_normal=difference_normal-alpha*vec(4)
  dd=0.0d0
  do i=1,3
   beta=inner_product_state(left_difference_normal, &
                               lleft_eigenvectors(i))
   beta=beta/inner_product_state(vec(i),lleft_eigenvectors(i))
   beta=dsign(1.0d0,beta)*dmin1(dabs(beta),limiter)
   dd=dd+beta*vec(i)
  end do
  left_difference_normal=dd    
  call transform_2(left_difference_normal,nml,left_difference)
  l_updated_state=l_updated_state+left_difference

 case(2)
! First, update the left and right states by adjusting the normal 
! velocity and pressure.
  ordin_state=p_cell%or_state
  call find_physical_state(l_updated_state,nml,left_phy_state)
  l_normal_velo=left_phy_state%value(2)
  l_pressure=left_phy_state%value(4)
  call find_physical_state(r_updated_state,nml,right_phy_state)
  r_normal_velo=right_phy_state%value(2)
  r_pressure=right_phy_state%value(4)
  call find_physical_state(ordin_state,nml,ordin_phy_state)
  or_normal_velo=ordin_phy_state%value(2)
  or_pressure=ordin_phy_state%value(4) 
  
  ttt=or_normal_velo-l_normal_velo
  ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),limiter)
  l_normal_velo=l_normal_velo+ttt
  ttt=or_normal_velo-r_normal_velo
  ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),limiter)
  r_normal_velo=r_normal_velo+ttt
   
  ttt=or_pressure-l_pressure
  ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),limiter)
  l_pressure=l_pressure+ttt
  ttt=or_pressure-r_pressure
  ttt=dsign(1.0d0,ttt)*dmin1(dabs(ttt),limiter)
  r_pressure=r_pressure+ttt
  
  left_phy_state%value(2)=l_normal_velo
  left_phy_state%value(4)=l_pressure
  right_phy_state%value(2)=r_normal_velo
  right_phy_state%value(4)=r_pressure
  
  call find_conservative_state(left_phy_state,nml,l_updated_state)
  call find_conservative_state(right_phy_state,nml,r_updated_state)
   
! The update the left and right states by enforcing the conservation
! relation.
  select case(g_cell%g_type)
   case('xx ')
   if(g_cell%point(1).eq.'left  ') then
    position=g_cell%mdis_posi(1)
   else
    position=-g_cell%mdis_posi(1)
   end if
   case('yy ')
    if(g_cell%point(1).eq.'left  ') then
	 position=g_cell%mdis_posi(2)
    else
     position=-g_cell%mdis_posi(2)
    end if	 	 
   case default; call error_message
  end select
  call transform_1(l_updated_state,nml,lnu)
  call transform_1(r_updated_state,nml,rnu)
  call transform_1(ordin_state,nml,dd)
  dd=dd-rnu-(position+0.5d0)*(lnu-rnu)
!  do i=1,4
!   lnu%value(i)=lnu%value(i)+dsign(1.0d0,dd%value(i))*dmin1(dabs(dd%value(i)),0.5d0*limiter) 
!   rnu%value(i)=rnu%value(i)+dsign(1.0d0,dd%value(i))*dmin1(dabs(dd%value(i)),0.5d0*limiter)
!  end do
  lnu%value(3)=lnu%value(3)+dsign(1.0d0,dd%value(3))*dmin1(dabs(dd%value(3)),0.1d0*limiter) 
  rnu%value(3)=rnu%value(3)+dsign(1.0d0,dd%value(3))*dmin1(dabs(dd%value(3)),0.1d0*limiter)
  call transform_2(lnu,nml,l_updated_state)
  call transform_2(rnu,nml,r_updated_state) 
  
end select 

c_cell%l_state(0)=l_updated_state
c_cell%r_state(0)=r_updated_state

! call transform_1(l_updated_state,nml,xxx)
! call transform_2(xxx,nml,yyy)
! continue

contains


subroutine form_intermediate

implicit none

real(8), dimension(4) :: wwl, wwr, ww
real(8) :: rhol, unl, utl, pl, rhor, unr, utr, pr, rho, un, ut, p  

wwl=left_nearest_in_normal%value
rhol=wwl(1)
unl=uf(wwl)
utl=vf(wwl)
pl=pf(wwl)
wwr=right_nearest_in_normal%value
rhor=wwr(1)
unr=uf(wwr)
utr=vf(wwr)
pr=pf(wwr)
rho=rhol
un=unr
ut=utr
p=pr
call tran3(rho,un,ut,p,ww)
intermediate_nearest_in_normal%value=ww

end subroutine form_intermediate


subroutine produce_limiter
! This subroutine produces the limiter.

implicit none

real(8), dimension(4) :: vv
real(8) :: aaa, bbb
integer :: i

aaa=0.0d0
do i=1,4
 vv(i)=((p_cell%l_state%value(i)-p_cell%r_state%value(i)))
 bbb=dabs(p_cell%l_state%value(i)+p_cell%r_state%value(i))
 if(bbb.gt.1.0d-6) then
  vv(i)=vv(i)/(p_cell%l_state%value(i)+p_cell%r_state%value(i))
 else
  vv(i)=0.0d0
 end if
 vv(i)=dabs(vv(i))
 aaa=aaa+vv(i)
end do
limiter=0.0d0
do i=1, 4
 limiter=limiter+vv(i)/aaa*dabs(difference%value(i))
end do
limiter=dsqrt(limiter)
limiter=limiter*limiter*limiter 

end subroutine


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
call roe(w2,w1,el,er,ch,coef)
do i=1,4
 left_eigenvectors(i)%value=el(i,:)
 right_eigenvectors(i)%value=er(i,:)
end do

end subroutine roe_decomposition


end subroutine evaluate_updated_state


subroutine produce_integral_points(side)

implicit none

character*6, intent(in) :: side
! The side the integral points are produced.

integer :: i, j
real(8), dimension(-1:1) :: xx, yy

select case(side)
 case('both  ')
  do i=-1,1
   do j=-1,1
    integral_points(i,j,1)=dfloat(i)*0.5d0
    integral_points(i,j,2)=dfloat(j)*0.5d0
   end do
  end do
  
 case('left  ','right ')
  
  select case(g_cell%g_type)
   case('xx ')
    if(g_cell%edge(1).eq.1) then
     xx=x
    else
     xx(-1)=x(1); xx(0)=x(0); xx(1)=x(-1)
    end if
    do j=-1, 1
     if(g_cell%point(1).eq.side) then
      integral_points(-1,j,1)=-0.5d0
      integral_points(1,j,1)=xx(j)
      integral_points(0,j,1)=0.5d0*(-0.5d0+xx(j))
     else
      integral_points(-1,j,1)=xx(j)
      integral_points(1,j,1)=0.5d0
      integral_points(0,j,1)=0.5d0*(xx(j)+0.5d0)
     end if
     do i=-1,1
      integral_points(i,j,2)=dfloat(j)*0.5d0
     end do
    end do
   case('yy ')
    if(g_cell%edge(1).eq.4) then
     yy=y
    else
     yy(-1)=y(1); yy(0)=y(0); yy(1)=y(-1)
    end if
	do i=-1,1
     if(g_cell%point(1).eq.side) then
      integral_points(i,-1,2)=-0.5d0
      integral_points(i,1,2)=yy(i)
      integral_points(i,0,2)=0.5d0*(-0.5d0+yy(i))
     else
      integral_points(i,-1,2)=yy(i)
      integral_points(i,1,2)=0.5d0
      integral_points(i,0,2)=0.5d0*(yy(i)+0.5d0)
     end if
     do j=-1,1
      integral_points(i,j,1)=dfloat(j)*0.5d0
     end do
    end do
  end select
    
end select 

end subroutine produce_integral_points


subroutine integral_state(side,original_or_updated,new_state)
! This subroutine evaluate the updated state on a given side 'side'.

implicit none
character*6, intent(in) :: side
character*8, intent(in) :: original_or_updated
type(state), intent(out) :: new_state

type(state), dimension(-1:1,0:2) :: inter_data
real(8), dimension(-1:1,0:2) :: inter_pts
integer :: k, i, j
logical :: ok

type(state), dimension(0:1) :: interpolation_data_1
real(8), dimension(0:1) :: interpolation_points_1
type(state), dimension(0:0) :: interpolation_data_2
real(8), dimension(0:0) :: interpolation_points_2
!type(state) :: constant_state

!constant_state=error_data
do k=-1,1
 call produce_interpolation_points(side,k)
 call produce_interpolation_data(side,k,original_or_updated)
 inter_data(k,:)=interpolation_data
 inter_pts(k,:)=interpolation_points
 
! if(k.eq.0) constant_state=interpolation_data(0)

 select case(g_cell%g_type)
  case('xx ')
   do i=-1,1
    integral_data(i,k)=lagrange_state(interpolation_data, &
                       interpolation_points,integral_points(i,k,1),2)
   end do
  case('yy ')
   do j=-1,1
    integral_data(k,j)=lagrange_state(interpolation_data, &
                       interpolation_points,integral_points(k,j,2),2)
   end do
 end select
end do

! The following checks whether the second-orderly computed 
! interpolation_data are OK.
do i=-1,1
 do j=-1,1
  ok=whether_state_OK(integral_data(i,j))
  if(.not.ok) exit
 end do
 if(.not.ok) exit
end do

if(.not.ok) then
 do k=-1,1
  interpolation_data_1=inter_data(k,0:1)
  interpolation_points_1=inter_pts(k,0:1)

  select case(g_cell%g_type)
   case('xx ')
    do i=-1,1
     integral_data(i,k)=lagrange_state(interpolation_data_1, &
                       interpolation_points_1,integral_points(i,k,1),1)
    end do
   case('yy ')
    do j=-1,1
     integral_data(k,j)=lagrange_state(interpolation_data_1, &
                        interpolation_points_1,integral_points(k,j,2),1)
    end do
  end select
 end do

! The following checks whether the first-orderly computed 
! interpolation_data are OK.
 do i=-1,1
  do j=-1,1
   ok=whether_state_OK(integral_data(i,j))
   if(.not.ok) exit
  end do
  if(.not.ok) exit 
 end do
end if

if(.not.ok) then
 do k=-1,1
  interpolation_data_2=inter_data(k,0:0)
  interpolation_points_2=inter_pts(k,0:0)

  select  case(g_cell%g_type)
   case('xx ')
    do i=-1,1
     integral_data(i,k)=lagrange_state(interpolation_data_2, &
                        interpolation_points_2,integral_points(i,k,1),0)
!     integral_data(i,k)=interpolation_data(0)
    end do
   case('yy ')
    do j=-1,1
     integral_data(k,j)=lagrange_state(interpolation_data_2, &
                        interpolation_points_2,integral_points(k,j,2),0)
!     integral_data(k,j)=interpolation_data(0)
    end do
  end select
 end do
end if

new_state=simpson_integral(g_cell%g_type)

end subroutine integral_state


function simpson_integral(fashion) result(c)

implicit none
character*3, intent(in) :: fashion
type(state) :: c

type(state) us1, us2, us3, simpson_state

select case(fashion)
 case('xx ')
  us1=integral_data(-1,-1)+4.0d0*integral_data(0,-1)+ &
      integral_data(1,-1)
  us1=((integral_points(1,-1,1)-integral_points(-1,-1,1))/6.0d0)*us1
  us2=integral_data(-1,0)+4.0d0*integral_data(0,0)+integral_data(1,0)
  us2=((integral_points(1,0,1)-integral_points(-1,0,1))/6.0d0)*us2
  us3=integral_data(-1,1)+4.0d0*integral_data(0,1)+integral_data(1,1)
  us3=((integral_points(1,1,1)-integral_points(-1,1,1))/6.0d0)*us3
 case('yy ')
  us1=integral_data(-1,-1)+4.0d0*integral_data(-1,0)+ &
      integral_data(-1,1)
  us1=((integral_points(-1,1,2)-integral_points(-1,-1,2))/6.0d0)*us1
  us2=integral_data(0,-1)+4.0d0*integral_data(0,0)+integral_data(0,1)
  us2=((integral_points(0,1,2)-integral_points(0,-1,2))/6.0d0)*us2
  us3=integral_data(1,-1)+4.0d0*integral_data(1,0)+integral_data(1,1)
  us3=((integral_points(1,1,2)-integral_points(1,-1,2))/6.0d0)*us3
end select
simpson_state=us1+4.0d0*us2+us3
simpson_state=simpson_state/6.0d0

c=simpson_state

end function simpson_integral


subroutine produce_interpolation_points(side,index)
! This subroutine produces the integral points on a give side 'side'.

implicit none
character*6, intent(in) :: side
integer, intent(in) :: index
real(8) :: step
integer :: i

if(g_cell%point(1).eq.side) then
 step=-0.50d0
else 
 step=0.5d0
end if
select case(g_cell%g_type)
 case('xx ')
!  interpolation_points(0)=x(index)
  if(g_cell%edge(1).eq.1) then
   interpolation_points(0)=x(index)
  else
   interpolation_points(0)=x(-index)
  end if
 case('yy ')
!  interpolation_points(0)=y(index)
  if(g_cell%edge(1).eq.4) then
   interpolation_points(0)=y(index)
  else
   interpolation_points(0)=y(-index)
  end if    
end select
do i=1,2
 interpolation_points(i)=interpolation_points(0)+ &
 (dfloat(i)+0.5d0)*step
end do

end subroutine produce_interpolation_points


subroutine produce_interpolation_data(side,index,original_or_updated)
! This subroutine produces the integral data on a give side 'side'.

implicit none
character*6, intent(in) :: side
integer, intent(in) :: index
character*8, intent(in) :: original_or_updated

integer :: kk

select case(g_cell%g_type)
 case('xx ')
  do kk=0,2
   if(g_cell%edge(1).eq.1) then
    interpolation_data(kk)=p_side(temp,interpolation_points(kk),y(index),side)
   else
    interpolation_data(kk)=p_side(temp,interpolation_points(kk),y(-index),side)
   end if
  end do
 case('yy ')
  do kk=0,2
   if(g_cell%edge(1).eq.4) then
    interpolation_data(kk)=p_side(temp,x(index),interpolation_points(kk),side)
   else
    interpolation_data(kk)=p_side(temp,x(-index),interpolation_points(kk),side)
   end if
  end do
 case default; call error_message
end select

if(original_or_updated.eq.'updated ') then
 select case(side)

  case('left  ')
   select case(g_cell%g_type)
    case('xx ')
     if(g_cell%edge(1).eq.1) then
      interpolation_data(0)=left_nearest(index)
     else
      interpolation_data(0)=left_nearest(-index)
     end if    
    case('yy ')
     if(g_cell%edge(1).eq.4) then
      interpolation_data(0)=left_nearest(index)
     else
      interpolation_data(0)=left_nearest(-index)
     end if    
    case default; call error_message
   end select
  case('right ')
   select case(g_cell%g_type)
    case('xx ')
     if(g_cell%edge(1).eq.1) then
      interpolation_data(0)=right_nearest(index)
     else
      interpolation_data(0)=right_nearest(-index)
     end if    
    case('yy ')
     if(g_cell%edge(1).eq.4) then
      interpolation_data(0)=right_nearest(index)
     else
      interpolation_data(0)=right_nearest(-index)
     end if
    case default; call error_message
   end select
  case default; call error_message
 end select

end if

end subroutine produce_interpolation_data


end module evaluate_updated_xxyy
