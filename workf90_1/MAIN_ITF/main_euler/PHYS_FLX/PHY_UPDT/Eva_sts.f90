module evaluate_updated_xxyy

use produce_polynomial_4_xx_yy
! 'proc_pln.f90'

use euler_functions
! 'euler_functions.f90'

type(auxi_crit_cell), pointer :: temp, boundary_end
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
		l_updated_state, r_updated_state, evaluate_updated_state_2


contains


subroutine update_xx_yy(acv)
! Do the update in 'xx'- and 'yy'-type critical cells, which is 
! done by extrapolation and numerical integral.

implicit none
type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch

integer :: i
real(8), dimension(2) :: pt, nml
real(8) :: velocity, rrr 
character*12 :: wave_type
type(state) :: fff

! type(geo_info) :: gggxxx

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
   
  select case(p_cell%wv_nb)
     
   case(1,3)
     
    do i=-1,1
     l_state=p_side(temp,x(i),y(i),3,'left  ')
     r_state=p_side(temp,x(i),y(i),3,'right ')
     nml=normal(i,:)
     call riemann(l_state,r_state,nml,p_cell%wv_nb,x(i),y(i), &
                  left_nearest(i),right_nearest(i),velocity,wave_type)
     fff=left_nearest(i)
	 fff=right_nearest(i)
    end do
     
! Evaluate the updated reconstruction solution in the critical cell and store the data in
! 'c_cell'.
    call evaluate_updated_state
    
   case(2) 
    
    call evaluate_updated_state_2

    fff=c_cell%l_state(0)
    
   case default; call error_message
   
  end select
            
 end if

 temp%c_cell=c_cell

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine update_xx_yy


subroutine evaluate_updated_state
! It evaluates the updated left and right states. It consists of 
! two steps: integration and clean-up.

implicit none

type(state) :: l_original_state, r_original_state
type(state) :: fff

integer :: ii, jj
type(state) :: left_phy_state, right_phy_state, ordin_state, ordin_phy_state
real(8) :: l_tangential_velo, r_tangential_velo, or_tangential_velo
real(8) :: l_x_velo, r_x_velo, or_x_velo, l_y_velo, r_y_velo, or_y_velo
real(8) :: l_pressure, r_pressure, or_pressure, ttt

type(state), dimension(-1:1,-1:1) :: left_intdata, right_intdata, l_intda, r_intda
type(state) :: left_physics, right_physics, physics
logical :: left_other_side, right_other_side
real(8) :: xxyy

! First, produce tentative updated states.

call produce_integral_points('both  ')
l_original_state=p_cell%l_state
r_original_state=p_cell%r_state
 
select case(p_cell%wv_nb)
 
 case(1,3)
  call integral_state('left  ','updated ',l_updated_state)
  call integral_state('right ','updated ',r_updated_state)
 
 case(2)
  call integral_state('left  ','updated ',l_updated_state)
  left_intdata=integral_data; l_intda=integral_data 
  call integral_state('right ','updated ',r_updated_state)
  right_intdata=integral_data; r_intda=integral_data
  
  select case(g_cell%g_type)
   
   case('xx ')
    do jj=-1,1
     do ii=-1,1
      if(g_cell%edge(1).eq.1) then
       xxyy=x(jj)
      else
       xxyy=x(-jj)    
      end if
      if(g_cell%point(1).eq.'left  ') then
       if(dfloat(ii)*0.5d0.gt.xxyy) then
        left_other_side=.true.; right_other_side=.false.
       else
        left_other_side=.false.; right_other_side=.true.
       end if
      else
       if(dfloat(ii)*0.5d0.gt.xxyy) then
        left_other_side=.false.; right_other_side=.true.
       else
        left_other_side=.true.; right_other_side=.false.
       end if
      end if
      call find_physical_state(left_intdata(ii,jj),(/1.0d0,0.0d0/),left_physics)
      call find_physical_state(right_intdata(ii,jj),(/1.0d0,0.0d0/),right_physics)
      if(left_other_side) then
       physics%value(1)=left_physics%value(1)
       physics%value(2)=right_physics%value(2)
       physics%value(3)=right_physics%value(3)
       physics%value(4)=right_physics%value(4)
       physics%gamma=right_physics%gamma
       call find_conservative_state(physics,(/1.0d0,0.0d0/),l_intda(ii,jj))
      end if 
      if(right_other_side) then
       physics%value(1)=right_physics%value(1)
       physics%value(2)=left_physics%value(2)
       physics%value(3)=left_physics%value(3)
       physics%value(4)=left_physics%value(4)
       physics%gamma=left_physics%gamma
       call find_conservative_state(physics,(/1.0d0,0.0d0/),r_intda(ii,jj))
      end if
     end do
    end do
   
   case('yy ')
    do ii=-1,1
     do jj=-1,1
      if(g_cell%edge(1).eq.4) then
       xxyy=y(ii)
      else
       xxyy=y(ii)    
      end if
      if(g_cell%point(1).eq.'left  ') then
       if(dfloat(jj)*0.5d0.gt.xxyy) then
        left_other_side=.true.; right_other_side=.false.
       else
        left_other_side=.false.; right_other_side=.true.
       end if
      else
       if(dfloat(jj)*0.5d0.gt.xxyy) then
        left_other_side=.false.; right_other_side=.true.
       else
        left_other_side=.true.; right_other_side=.false.
       end if
      end if
      call find_physical_state(left_intdata(ii,jj),(/1.0d0,0.0d0/),left_physics)
      call find_physical_state(right_intdata(ii,jj),(/1.0d0,0.0d0/),right_physics)
      if(left_other_side) then
       physics%value(1)=left_physics%value(1)
       physics%value(2)=right_physics%value(2)
       physics%value(3)=right_physics%value(3)
       physics%value(4)=right_physics%value(4)
	   physics%gamma=right_physics%gamma
       call find_conservative_state(physics,(/1.0d0,0.0d0/),l_intda(ii,jj))
      end if 
      if(right_other_side) then
       physics%value(1)=right_physics%value(1)
       physics%value(2)=left_physics%value(2)
       physics%value(3)=left_physics%value(3)
       physics%value(4)=left_physics%value(4)
	   physics%gamma=right_physics%gamma
       call find_conservative_state(physics,(/1.0d0,0.0d0/),r_intda(ii,jj))
      end if
     end do
    end do
    
   case default; call error_message
    
  end select

!   do ii=-1,1
!    do jj=-1,1
!      call find_physical_state(l_intda(ii,jj),(/1.0d0,0.0d0/),left_physics)
!       call find_physical_state(r_intda(ii,jj),(/1.0d0,0.0d0/),right_physics)
!    end do
!   end do

  integral_data=l_intda
  l_updated_state=simpson_integral(g_cell%g_type)
  integral_data=r_intda
  r_updated_state=simpson_integral(g_cell%g_type)

!   call find_physical_state(l_updated_state,(/1.0d0,0.0d0/),left_physics)
!   call find_physical_state(r_updated_state,(/1.0d0,0.0d0/),right_physics)
!   call find_physical_state(p_cell%or_state,(/1.0d0,0.0d0/),physics)

 case default; call error_message
    
end select
  
c_cell%l_state(0)=l_updated_state
fff=c_cell%l_state(0)
c_cell%r_state(0)=r_updated_state
fff=c_cell%r_state(0)

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
type(state) :: fff

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
	fff=integral_data(i,k)
   end do
  case('yy ')
   do j=-1,1
    integral_data(k,j)=lagrange_state(interpolation_data, &
                       interpolation_points,integral_points(k,j,2),2)
    fff=integral_data(i,k)
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

type(state) :: us1, us2, us3, simpson_state

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
 step=-1.00d0
else 
 step=1.0d0
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
 interpolation_points(i)=interpolation_points(0)+dfloat(i)*step
end do

end subroutine produce_interpolation_points


subroutine produce_interpolation_data(side,index,original_or_updated)
! This subroutine produces the integral data on a give side 'side'.

implicit none
character*6, intent(in) :: side
integer, intent(in) :: index
character*8, intent(in) :: original_or_updated

integer :: kk
type(state), dimension(0:2) :: fff

select case(g_cell%g_type)
 case('xx ')
  do kk=0,2
   if(g_cell%edge(1).eq.1) then
    interpolation_data(kk)=p_side(temp,interpolation_points(kk), &
	                              y(index),3,side)
    fff=interpolation_data(kk)
   else
    interpolation_data(kk)=p_side(temp,interpolation_points(kk), &
	                              y(-index),3,side)
    fff=interpolation_data(kk)
   end if
  end do
 case('yy ')
  do kk=0,2
   if(g_cell%edge(1).eq.4) then
    interpolation_data(kk)=p_side(temp,x(index), &
	                              interpolation_points(kk),3,side)
	fff=interpolation_data(0)
   else
    interpolation_data(kk)=p_side(temp,x(-index), &  
	                              interpolation_points(kk),3,side)
	fff=interpolation_data(0)
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
      fff=interpolation_data(0)
     else
      interpolation_data(0)=left_nearest(-index)
      fff=interpolation_data(0)
     end if    
    case('yy ')
     if(g_cell%edge(1).eq.4) then
      interpolation_data(0)=left_nearest(index)
	  fff=interpolation_data(0)
     else
      interpolation_data(0)=left_nearest(-index)
	  fff=interpolation_data(0)
     end if    
    case default; call error_message
   end select
  case('right ')
   select case(g_cell%g_type)
    case('xx ')
     if(g_cell%edge(1).eq.1) then
      interpolation_data(0)=right_nearest(index)
	  fff=interpolation_data(0)
     else
      interpolation_data(0)=right_nearest(-index)
	  fff=interpolation_data(0)
     end if    
    case('yy ')
     if(g_cell%edge(1).eq.4) then
      interpolation_data(0)=right_nearest(index)
	  fff=interpolation_data(0)
     else
      interpolation_data(0)=right_nearest(-index)
	  fff=interpolation_data(0)
     end if
    case default; call error_message
   end select
  case default; call error_message
 end select

end if

end subroutine produce_interpolation_data


subroutine evaluate_updated_state_2
! It evaluates the updated left and right states for the second-field. 
! It consists of two steps: integration and clean-up.

implicit none

type(state), dimension(-1:1,-1:1) :: left_intdata, right_intdata, &
                                     l_intda, r_intda
type(state), dimension(-1:1,-1:1) :: left_physics, right_physics, &
                                     new_l_phy, new_r_phy
real(8), dimension(2) :: point
integer :: ii, jj
character*6 :: side
logical :: lower_order
type(state) :: fff

call produce_integral_points('both  ')

!Firstly, find the data on the left side.
do ii=-1,1
 do jj=-1,1
  left_intdata(ii,jj)=p_side(temp,integral_points(ii,jj,1), &
                             integral_points(ii,jj,2),3,'left  ')
  
  fff=left_intdata(ii,jj)

 end do
end do
do ii=-1,1
 do jj=-1,1
  call find_physical_state(left_intdata(ii,jj),(/1.0d0,0.0d0/), &
                           left_physics(ii,jj))
 
  fff=left_physics(ii,jj)
 
 end do
end do
lower_order=.false.
do ii=-1,1
 do jj=-1,1
  if(left_physics(ii,jj)%value(1).lt.0.0d0.or.left_physics(ii,jj)% &
     value(4).lt.0.0d0) then
   lower_order=.true.; exit
  end if
 end do
 if(lower_order) exit
end do

if(.not.lower_order) goto 10

do ii=-1,1
 do jj=-1,1
  left_intdata(ii,jj)=p_side(temp,integral_points(ii,jj,1), &
                             integral_points(ii,jj,2),2,'left  ')

  fff=left_intdata(ii,jj)

 end do
end do
do ii=-1,1
 do jj=-1,1
  call find_physical_state(left_intdata(ii,jj),(/1.0d0,0.0d0/), &
                           left_physics(ii,jj))

  fff=left_physics(ii,jj)

 end do
end do
lower_order=.false.
do ii=-1,1
 do jj=-1,1
  if(left_physics(ii,jj)%value(1).lt.0.0d0.or.left_physics(ii,jj)% &
     value(4).lt.0.0d0) then
   lower_order=.true.; exit
  end if
 end do
end do

if(.not.lower_order) goto 10

do ii=-1,1
 do jj=-1,1
  left_intdata(ii,jj)=p_side(temp,integral_points(ii,jj,1), &
                             integral_points(ii,jj,2),1,'left  ')

  fff=left_intdata(ii,jj)

 end do
end do
do ii=-1,1
 do jj=-1,1
  call find_physical_state(left_intdata(ii,jj),(/1.0d0,0.0d0/), &
                           left_physics(ii,jj))

  fff=left_physics(ii,jj)

 end do
end do

10 continue 

! Then, find the data on the right side.
do ii=-1,1
 do jj=-1,1
  right_intdata(ii,jj)=p_side(temp,integral_points(ii,jj,1), &
                            integral_points(ii,jj,2),3,'right ')
 
  fff=right_intdata(ii,jj)
 
 end do
end do
do ii=-1,1
 do jj=-1,1
  call find_physical_state(right_intdata(ii,jj),(/1.0d0,0.0d0/), &
                           right_physics(ii,jj))

  fff=right_physics(ii,jj)

 end do
end do
lower_order=.false.
do ii=-1,1
 do jj=-1,1
  if(right_physics(ii,jj)%value(1).lt.0.0d0.or.right_physics(ii,jj)% &
     value(4).lt.0.0d0) then
   lower_order=.true.; exit
  end if  
 end do
end do

if(.not.lower_order) goto 20

do ii=-1,1
 do jj=-1,1
  right_intdata(ii,jj)=p_side(temp,integral_points(ii,jj,1), &
                            integral_points(ii,jj,2),2,'right ')

  fff=right_intdata(ii,jj)

 end do
end do
do ii=-1,1
 do jj=-1,1
  call find_physical_state(right_intdata(ii,jj),(/1.0d0,0.0d0/), &
                           right_physics(ii,jj))

  fff=right_physics(ii,jj)

 end do
end do
lower_order=.false.
do ii=-1,1
 do jj=-1,1
  if(right_physics(ii,jj)%value(1).lt.0.0d0.or.right_physics(ii,jj)% &
     value(4).lt.0.0d0) then
   lower_order=.true.; exit
  end if  
 end do
end do

if(.not.lower_order) goto 20

do ii=-1,1
 do jj=-1,1
  right_intdata(ii,jj)=p_side(temp,integral_points(ii,jj,1), &
                            integral_points(ii,jj,2),1,'right ')

  fff=right_intdata(ii,jj)

 end do
end do
do ii=-1,1
 do jj=-1,1
  call find_physical_state(right_intdata(ii,jj),(/1.0d0,0.0d0/), &
                           right_physics(ii,jj))

  fff=right_physics(ii,jj)

 end do
end do

20 continue

! Finally, compose the data on the two side in keeping the velocity and
! pressure continuous.
new_l_phy%value(1)=left_physics%value(1)
new_r_phy%value(1)=right_physics%value(1)
do ii=-1,1
 do jj=-1,1
  point=integral_points(ii,jj,:)
  side=side_of_point(g_cell,point)  
  select case(side)
   case('left  ')
    new_l_phy(ii,jj)%value(2)=left_physics(ii,jj)%value(2)
    new_l_phy(ii,jj)%value(3)=left_physics(ii,jj)%value(3)
    new_l_phy(ii,jj)%value(4)=left_physics(ii,jj)%value(4)
	new_l_phy(ii,jj)%gamma=left_physics(ii,jj)%gamma
	fff=new_l_phy(ii,jj)
   case('right ')
    new_l_phy(ii,jj)%value(2)=right_physics(ii,jj)%value(2)
    new_l_phy(ii,jj)%value(3)=right_physics(ii,jj)%value(3)
    new_l_phy(ii,jj)%value(4)=right_physics(ii,jj)%value(4)
	new_l_phy(ii,jj)%gamma=left_physics(ii,jj)%gamma
! very important sentence

	fff=new_l_phy(ii,jj)
   case default; call error_message
  end select
  new_r_phy(ii,jj)%value(2)=new_l_phy(ii,jj)%value(2) 
  new_r_phy(ii,jj)%value(3)=new_l_phy(ii,jj)%value(3) 
  new_r_phy(ii,jj)%value(4)=new_l_phy(ii,jj)%value(4) 
  new_r_phy(ii,jj)%gamma=right_physics(ii,jj)%gamma
! very important sentence

  fff=new_r_phy(ii,jj)
 end do
end do

do ii=-1,1
 do jj=-1,1
  call find_conservative_state(new_l_phy(ii,jj),(/1.0d0,0.0d0/), &
                               left_intdata(ii,jj))

  fff=left_intdata(ii,jj)

  call find_conservative_state(new_r_phy(ii,jj),(/1.0d0,0.0d0/), &
                               right_intdata(ii,jj))

  fff=right_intdata(ii,jj)

 end do
end do

integral_data=left_intdata
l_updated_state=simpson_integral(g_cell%g_type)
integral_data=right_intdata
r_updated_state=simpson_integral(g_cell%g_type)

c_cell%l_state(0)=l_updated_state
fff=c_cell%l_state(0)
c_cell%r_state(0)=r_updated_state
fff=c_cell%r_state(0)

end subroutine evaluate_updated_state_2


end module evaluate_updated_xxyy
