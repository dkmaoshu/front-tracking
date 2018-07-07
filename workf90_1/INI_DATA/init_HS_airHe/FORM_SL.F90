module form_solution
! This module contains all the operations used in forming the 
! initial data of the solution.

use solution
! 'solution.f90'

use out_intermediate
! 'outputs.f90'

use grid_map
! 'grid_map.f90'

use numerical_integral_in_2D
! 'num_int2.f90'

use euler_functions
! 'euler_functions.f90'

implicit none

type(critical_cell), pointer :: current, temp
type(geo_info) :: g_cell, g_cell_1
type(phy_info) :: p_cell, p_cell_1
type(state) :: wt, wm, wb

real(8), dimension(3) :: values, points, coefficients
real(8), dimension(3,2) :: stencil
integer :: ii, jj, i, j, k, nn, first_idx, cases
real(8) :: pi
logical :: if_on_curve
real(8) :: xx, yy, dd, radius

public  set_in_smooth, set_discontinuity_curve, &
        set_slops_on_curve, compute_or_states, cases
private	f_1, f_2, f_3


contains


subroutine set_in_smooth

implicit none

real(8) :: rhot, ut, vt, pt, rhom, um, vm, pm, rhob, gammam, gammab, gammat
real(8), dimension(4) :: wwt, wwm, wwb

logical :: if_on_curve
real(8) :: ct, cb, cm, mu, m0, c_ast, dd
real(8), dimension(4) :: gt, gm
real(8) :: v_0, v_1, shock_speed, hc_check
integer :: i, j
type(state) :: ww


 rhom=1.225d0;  um=0.0d0;       vm=0.0d0; pm=1.01325d0; gammam=1.4d0
! rhot=1.6861d0; ut=-0.494137d0; vt=0.0d0; pt=2.50638d0; gammat=1.4d0
 rhob=0.2228d0;                                         gammab=1.648d0
 gammat=gammam


! The constant $\mu$.
  mu=dsqrt((gammam-1.0d0)/(gammam+1.0d0))


! The shock on the top is of Mach 5.
  m0=1.22

  cm=dsqrt(gammam*pm/rhom)

! $v_0$ is the relative velocity of the flow ahead of the shock (middle region).
  v_0=m0*cm

! $c_{\ast}$ is the critical speed of the flow across the shock.
  c_ast=mu*mu*v_0*v_0+(1.0d0-mu*mu)*cm*cm
  c_ast=dsqrt(c_ast)

  pt=((1.0d0+mu*mu)*m0*m0-mu*mu)*pm

! $v_1$ is the relative velocity of the flow behind the shock (top region).
  v_1=c_ast*c_ast/v_0

  ct=(c_ast*c_ast-mu*mu*v_1*v_1)/(1.0d0-mu*mu)
  ct=dsqrt(ct)

  rhot=gammam*pt/(ct*ct)

  shock_speed=-v_0
  ut=v_1-v_0

  um=um+0.4d0
  ut=ut+0.4d0


 call tran3(rhom,um,vm,pm,wwm,gammam)
 call tran3(rhot,ut,vt,pt,wwt,gammam)
 call tran3(rhob,um,vm,pm,wwb,gammab)

! wwm(1)=1.0d0
! wwm(2)=0.0d0
! wwm(3)=0.0d0
! wwm(4)=2.5d1
! wwt(1)=1.0d0
! wwt(2)=0.0d0
! wwt(3)=0.0d0
! wwt(4)=2.5d0

! gammam=1.09d0
! gammat=1.09d0
 
 wm%value=wwm
 wm%gamma=gammam

 wt%value=wwt
 wt%gamma=gammam

 wb%value=wwb
 wb%gamma=gammab

do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  call check_on_curve(i,j,if_on_curve)
  if(.not.if_on_curve) then
   dd=dsqrt((dfloat(i)*h)**2.0d0+(dfloat(j)*h)**2.0d0)
   xx=dfloat(i)*h; yy=dfloat(j)*h
   if(dd.gt.radius) then
     if(xx.lt.0.3d0-0.5d0*h) then
	  uu(i,j)%value=wm%value
	  uu(i,j)%gamma=gammam
     else
	  if(xx.gt.0.3d0+0.5d0*h) then
       uu(i,j)%value=wwt
	   uu(i,j)%gamma=gammat
      else	  
	   uu(i,j)%value=0.5d0*(wwm+wwt)
	   uu(i,j)%gamma=gammam
	  end if
	 end if  
   else
    uu(i,j)%value=wwb
    uu(i,j)%gamma=gammab 
   end if 
  end if
 end do
end do


contains


subroutine check_on_curve(x_idx,y_idx,if_on_curve)

implicit none
integer, intent(in) :: x_idx, y_idx
logical, intent(out) :: if_on_curve

real(8) :: x, y, rr, ss

if_on_curve=.false.
x=dfloat(x_idx)*h; y=dfloat(y_idx)*h
rr=dsqrt(x*x+y*y)
ss=0.5d0*dsqrt(2.0d0)*h
if(rr.ge.radius+ss) return
if(rr.le.radius-ss) return

if(dabs(y).gt.dabs(x)) then
 xx=x-0.5d0*h
 if(dabs(xx).le.radius) then
  yy=dsqrt(radius**2.0d0-xx**2.0d0)
  yy=(yy-dabs(y))/h
  if(dabs(yy).lt.0.5d0) then
   if_on_curve=.true.
   return
  end if
 end if
 xx=x+0.5d0*h
 if(dabs(xx).le.radius) then
  yy=dsqrt(radius**2.0d0-xx**2.0d0)
  yy=(yy-dabs(y))/h
  if(dabs(yy).lt.0.5d0) then
   if_on_curve=.true.
   return
  end if
 end if
else
 yy=y-0.5d0*h
 if(dabs(yy).le.radius) then
  xx=dsqrt(radius**2.0d0-yy**2.0d0)
  xx=(xx-dabs(x))/h
  if(dabs(xx).lt.0.5d0) then
   if_on_curve=.true.
   return
  end if
 end if
 yy=y+0.5d0*h
 if(dabs(yy).le.radius) then
  xx=dsqrt(radius**2.0d0-yy**2.0d0)
  xx=(xx-dabs(x))/h
  if(dabs(xx).lt.0.5d0) then
   if_on_curve=.true.
   return
  end if
 end if
end if

end subroutine check_on_curve


end subroutine set_in_smooth


subroutine set_discontinuity_curve(curve_number)

implicit none
integer, intent(in) :: curve_number

type(state) :: l_state, r_state, c1, c2
real(8), dimension(2) :: normal, pt1, pt2
real(8) :: velocity
character*12 :: wave_type

real(8) :: xposi, yposi

if(curve_number.eq.1) call clean_up_grid_map(ggd_cell)
print*, 'Creat the',curve_number,'th  discontinuity curve.'

cvv(curve_number)%status='awake'
cvv(curve_number)%cv_type='circular'
cvv(curve_number)%begin_end=0; cvv(curve_number)%end_end=0
allocate(cvv(curve_number)%begin)
nullifY(cvv(curve_number)%begin%previous)
current=>cvv(curve_number)%begin

! The first critical cell is set ot be on the (positive) x-axis.
ii=idint(radius/h+0.5d0); jj=0; first_idx=ii
g_cell%x_idx=ii; g_cell%y_idx=jj
yy=-0.5d0*h
xx=dsqrt(radius*radius-yy*yy)
g_cell%edge(1)=1
xx=xx/h-dfloat(g_cell%x_idx)
g_cell%dis_posi(1)=xx
yy=0.5d0*h
xx=dsqrt(radius*radius-yy*yy)
xx=xx/h-dfloat(g_cell%x_idx)
g_cell%edge(2)=3; g_cell%dis_posi(2)=xx
call find_points_from_edges(g_cell)
call find_type_from_edges(g_cell)

call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2)
normal=normal_of_line(pt1,pt2)
xx=dfloat(g_cell%x_idx)*h; yy=dfloat(g_cell%y_idx)*h
select case(cases)
 case(1); l_state=f_2(xx,yy)
 case(2); l_state=f_3(xx,yy,radius)
 case default; call error_message
end select
r_state=f_1(xx,yy)

xposi=dfloat(ii); yposi=dfloat(jj)
call riemann(l_state,r_state,normal,2,xposi,yposi,c1,c2,velocity,wave_type)
p_cell%l_state=c1; p_cell%r_state=c2

if(curve_number.eq.1) p_cell%l_state=l_state
if(curve_number.eq.1) p_cell%r_state=r_state


p_cell%wv_nb=curve_number+1

k=1
allocate(temp)
temp%g_cell=g_cell; temp%p_cell=p_cell
temp%address%idx=k; k=k+1
temp%address%cv_nmb=curve_number
select case(curve_number)
 case(1)
  temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)
 case(2)
  temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)
end select

uu(ii,jj)=error_data
ggd_cell(ii,jj)%region='crit'
do j=1,4
! select case(curve_number)
!  case(1)
   if(g_cell%point(j).eq.'left  ') then
    ggd_cell(ii,jj)%ccpt(j)%address=temp%address
    ggd_cell(ii,jj)%ccpt(j)%side=g_cell%point(j)
   end if
!  case(2)
   if(g_cell%point(j).eq.'right ') then
    ggd_cell(ii,jj)%ccpt(j)%address=temp%address
    ggd_cell(ii,jj)%ccpt(j)%side=g_cell%point(j)
   end if
! end select
end do

temp%previous=>current
current%next=>temp
current=>temp

do while(ii.ne.first_idx.or.jj.ne.-1)
 select case(g_cell%edge(2))
  case(1); jj=jj-1
  case(2); ii=ii+1
  case(3); jj=jj+1
  case(4); ii=ii-1
 end select 
 g_cell%x_idx=ii; g_cell%y_idx=jj
 g_cell%edge(1)=neighbor_edge(g_cell%edge(2))
 g_cell%dis_posi(1)=g_cell%dis_posi(2)

 do i=1,4
  if(i.eq.g_cell%edge(1)) cycle
  select case(i)
   case(1,3)
    if(i.eq.1) then
     yy=(dfloat(jj)-0.5d0)*h
    else
     yy=(dfloat(jj)+0.5d0)*h
    end if
    if(dabs(yy).gt.radius) cycle
    xx=dsqrt(radius*radius-yy*yy)
    if(ii.lt.0) then
     xx=-xx
    else
     if(ii.eq.0.and.jj.lt.0) xx=-xx
    end if
    xx=xx/h-dfloat(g_cell%x_idx)
    if(dabs(xx).lt.0.5d0) exit
   case(2,4)
    if(i.eq.2) then
     xx=(dfloat(ii)+0.5d0)*h
    else
     xx=(dfloat(ii)-0.5d0)*h
    end if
    if(dabs(xx).gt.radius) cycle
    yy=dsqrt(radius*radius-xx*xx)
    if(jj.lt.0) then
     yy=-yy
    else
     if(jj.eq.0.and.ii.lt.0) yy=-yy
    end if
	yy=yy/h-dfloat(g_cell%y_idx)
    if(dabs(yy).lt.0.5d0) exit
  end select
 end do
 if(i.gt.4) call error_message
 g_cell%edge(2)=i
 select case(i)
  case(1,3); g_cell%dis_posi(2)=xx
  case(2,4); g_cell%dis_posi(2)=yy
 end select     
 call find_points_from_edges(g_cell)
 call find_type_from_edges(g_cell)

 call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2)
 normal=normal_of_line(pt1,pt2)
 xx=dfloat(g_cell%x_idx)*h; yy=dfloat(g_cell%y_idx)*h
 select case(cases)
  case(1); l_state=f_2(xx,yy)
  case(2); l_state=f_3(xx,yy,radius)
  case default; call error_message
 end select
 r_state=f_1(xx,yy)
 call riemann(l_state,r_state,normal,2,dfloat(ii)*h,dfloat(jj)*h,c1,c2,velocity,wave_type)
 p_cell%l_state=c1; p_cell%r_state=c2
 
 if(curve_number.eq.1) p_cell%l_state=l_state
 if(curve_number.eq.1) p_cell%r_state=r_state
 

 p_cell%wv_nb=curve_number+1

 allocate(temp)
 temp%g_cell=g_cell; temp%p_cell=p_cell
 temp%address%idx=k; k=k+1
 temp%address%cv_nmb=curve_number
 select case(curve_number)
  case(1)
   temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)
  case(2)
   temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)
 end select

 uu(ii,jj)=error_data
 ggd_cell(ii,jj)%region='crit'
 do j=1,4
!  select case(curve_number)
!   case(1)
    if(g_cell%point(j).eq.'left  ') then
     ggd_cell(ii,jj)%ccpt(j)%address=temp%address
     ggd_cell(ii,jj)%ccpt(j)%side=g_cell%point(j)
    end if
!   case(2)
    if(g_cell%point(j).eq.'right ') then
     ggd_cell(ii,jj)%ccpt(j)%address=temp%address
     ggd_cell(ii,jj)%ccpt(j)%side=g_cell%point(j)
    end if
!  end select
 end do

 temp%previous=>current
 current%next=>temp
 current=>temp
end do

nullify(cvv(curve_number)%begin%next%previous)
allocate(cvv(curve_number)%eend)	
nullify(cvv(curve_number)%eend%next)
cvv(curve_number)%eend%previous=>current
cvv(curve_number)%total=k-1
cvv(curve_number)%wave=curve_number+1

cvv(curve_number)%begin%next%previous=>cvv(curve_number)%eend%previous
cvv(curve_number)%eend%previous%next=>cvv(curve_number)%begin%next

print*, 'The',curve_number,'th discontinuity curve has been created.' 
! Created the 'curve_number'th discontinuity curve.

end subroutine set_discontinuity_curve


subroutine set_slops_on_curve(curve_number)

implicit none
integer, intent(in) :: curve_number

integer :: head_mark, io, jo
logical :: head_switch
type(state), dimension(-3:3) :: uu

! type(state) :: xxx

nullify(temp)
if(associated(cvv(curve_number)%begin%next)) then
 temp=>cvv(curve_number)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch) 
 g_cell=temp%g_cell; p_cell=temp%p_cell
 temp%udiff_xl=error_data; temp%udiff_yl=error_data
 temp%udiff_xr=error_data; temp%udiff_yr=error_data

! The following computes the difference quotients in $x$- and $y$-directions.
 call smth_dpr(io,jo,temp,uu,'xx ','left  ','critical',.true.,2,'full    ')
 temp%udiff_xl=0.5d0*(uu(1)-uu(-1)) !/h

!  xxx=temp%udiff_xl

 call smth_dpr(io,jo,temp,uu,'yy ','left  ','critical',.true.,2,'full    ')
 temp%udiff_yl=0.5d0*(uu(1)-uu(-1)) !/h
 call smth_dpr(io,jo,temp,uu,'xx ','right ','critical',.true.,2,'full    ')
 temp%udiff_xr=0.5d0*(uu(1)-uu(-1)) !/h
 call smth_dpr(io,jo,temp,uu,'yy ','right ','critical',.true.,2,'full    ')
 temp%udiff_yr=0.5d0*(uu(1)-uu(-1)) !/h

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if

end do

end subroutine set_slops_on_curve


subroutine compute_or_states(curve_number)

implicit none
integer, intent(in) :: curve_number

integer :: head_mark
logical :: head_switch
type(geo_info) :: gp_cell, gn_cell
real(8), dimension(2) :: pt

nullify(temp)
if(associated(cvv(curve_number)%begin%next)) then
 temp=>cvv(curve_number)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

do while(associated(temp).and.head_switch) 
 g_cell=temp%g_cell; p_cell=temp%p_cell
 gp_cell=temp%previous%g_cell; gn_cell=temp%next%g_cell

! First form the stencil
 call pick_pt(g_cell,pt)
 stencil(2,:)=pt
 call pick_pt(gp_cell,pt)
 pt=pt+(/dfloat(gp_cell%x_idx-g_cell%x_idx), &
        dfloat(gp_cell%y_idx-g_cell%y_idx)/)
 stencil(1,:)=pt
 call pick_pt(gn_cell,pt)
 pt=pt+(/dfloat(gn_cell%x_idx-g_cell%x_idx), &
        dfloat(gn_cell%y_idx-g_cell%y_idx)/)
 stencil(3,:)=pt

! The following computes the ordinary cell-average. First, 
! determine the 2D integrant.
 u_constant=p_cell%r_state-p_cell%l_state
 slop_x=temp%udiff_xr-temp%udiff_xl
 slop_y=temp%udiff_yr-temp%udiff_yl

 a_upper=error_data; b_upper=error_data; c_upper=error_data
 a_lower=error_data; b_upper=error_data; c_lower=error_data

 select case(g_cell%g_type)
  case('xx ','yy '); call case_xx_yy
  case('xy ');       call case_xy
  case default; call error_message
 end select

 p_cell%or_state%gamma=error_data
 
 temp%p_cell=p_cell
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if

end do


contains


subroutine case_xx_yy

implicit none

type(state) :: difference

geometrical_type=g_cell%g_type

select case(geometrical_type)
 case('xx ')
  values=stencil(:,1); points=stencil(:,2)
 case('yy ')
  values=stencil(:,2); points=stencil(:,1)
end select
call polynomial(values,points,coefficients)

a_upper=coefficients(3); b_upper=coefficients(2)
c_upper=coefficients(1)
a_lower=0.0d0; b_lower=0.0d0; c_lower=-0.5d0
lower=-0.5d0; upper=0.5d0
call integral_2(difference)

select case(g_cell%point(1))
 case('left  '); p_cell%or_state=p_cell%r_state-difference
 case('right '); p_cell%or_state=p_cell%l_state+difference
 case default; call error_message
end select


end subroutine case_xx_yy


subroutine case_xy

implicit none

type(state) :: difference
integer :: single, edge_num
real(8), dimension(2) :: pt1, pt2, pt_corner
real(8) :: length1, length2

! Determine the geometrical type
call find_single(g_cell,single)
call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2)
call pick_corner(pt_corner,single)
length1=length_of_segment(pt1,pt_corner)
length2=length_of_segment(pt2,pt_corner)
if(length1.gt.length2) then
 edge_num=1
else
 edge_num=2
end if
select case(g_cell%edge(edge_num))
 case(1,3); geometrical_type='yy '
 case(2,4); geometrical_type='xx '
 case default; call error_message
end select

! Determine the curve segment.
select case(geometrical_type)
 case('xx ')
  values=stencil(:,1); points=stencil(:,2)
 case('yy ')
  values=stencil(:,2); points=stencil(:,1)
end select
call polynomial(values,points,coefficients)

! Determine integral limits.
select case(single)
 case(1)
  a_upper=coefficients(3); b_upper=coefficients(2)
  c_upper=coefficients(1)
  a_lower=0.0d0; b_lower=0.0d0; c_lower=-0.5d0
  lower=-0.5d0; upper=g_cell%dis_posi(edge_num)
 case(2)
  select case(geometrical_type)
   case('xx ')
    a_lower=coefficients(3); b_lower=coefficients(2)
    c_lower=coefficients(1)
    a_upper=0.0d0; b_upper=0.0d0; c_upper=0.5d0
    lower=-0.5d0; upper=g_cell%dis_posi(edge_num)
   case('yy ')
    a_upper=coefficients(3); b_upper=coefficients(2)
    c_upper=coefficients(1)
    a_lower=0.0d0; b_lower=0.0d0; c_lower=-0.5d0
    lower=g_cell%dis_posi(edge_num); upper=0.5d0
  end select
 case(3)
  a_lower=coefficients(3); b_lower=coefficients(2)
  c_lower=coefficients(1)
  a_upper=0.0d0; b_upper=0.0d0; c_upper=0.5d0
  lower=g_cell%dis_posi(edge_num); upper=0.5d0
 case(4)
  select case(geometrical_type)
   case('xx ')
    a_upper=coefficients(3); b_upper=coefficients(2)
    c_upper=coefficients(1)
    a_lower=0.0d0; b_lower=0.0d0; c_lower=-0.5d0
    lower=g_cell%dis_posi(edge_num); upper=0.5d0
   case('yy ')
    a_lower=coefficients(3); b_lower=coefficients(2)
    c_lower=coefficients(1)
    a_upper=0.0d0; b_upper=0.0d0; c_upper=0.5d0
    lower=-0.5d0; upper=g_cell%dis_posi(edge_num)
  end select
end select

! Compute the difference.
call integral_2(difference)

select case(g_cell%point(single))
 case('left  '); p_cell%or_state=p_cell%r_state-difference
 case('right '); p_cell%or_state=p_cell%l_state+difference
 case default; call error_message
end select

end subroutine case_xy


end subroutine compute_or_states


function f_1(x,y) result(c)

implicit none
real(8), intent(in) :: x, y
type(state) :: c

real(8) :: rhot, ut, vt, pt, rhom, um, vm, pm, rhob, gammam, gammab, gammat
real(8), dimension(4) :: wwt, wwm, wwb


real(8) :: a
a=x+y

rhom=1.225d0;  um=0.4d0;       vm=0.0d0; pm=1.01325d0; gammam=1.4d0
call tran3(rhom,um,vm,pm,wwm,gammam)

!c%value(1)=1.0d0; c%value(2)=0.0d0
!c%value(3)=0.0d0; c%value(4)=1.0d0/9.0d-2
c%value=wwm
c%gamma=gammam

end function f_1


function f_2(x,y) result(c)

implicit none
real(8), intent(in) :: x, y
type(state) :: c
real(8) :: rhot, ut, vt, pt, rhom, um, vm, pm, rhob, gammam, gammab, gammat
real(8), dimension(4) :: wwt, wwm, wwb

!real(8) :: pi
real(8) :: a

a=x+y
!pi=4.0d0*datan(1.0d0)
!c=2.0d0+0.75d0*dcos((x+y)*(x+y)*pi)*dcos((x-y)*(x-y)*pi)
!c=2.0d0+0.75d0*dcos((x+y)*pi)*dcos((x-y)*pi)

rhob=0.2228d0;  um=0.4d0; vm=0.0d0; pm=1.01325d0;  gammab=1.648d0
call tran3(rhob,um,vm,pm,wwb,gammab)

!c%value(1)=1.2d0; c%value(2)=0.0d0
!c%value(3)=0.0d0; c%value(4)=1.0d0/9.0d-2
c%value=wwb
c%gamma=gammab

end function f_2


function f_3(x,y,radius) result(c)

implicit none
real(8), intent(in) :: x, y, radius
type(state) :: c

real(8) :: rhot, ut, vt, pt, rhom, um, vm, pm, rhob, gammam, gammab, gammat
real(8), dimension(4) :: wwt, wwm, wwb


real(8) :: rr, uu, vv, rhoo, pp, cos, sin, uuu, gg
real(8), dimension(4) :: ww

!rhoo=0.62847d0; uuu=1.6596d0; pp=5.2191d0
!gg=1.09d0

!if(dabs(x).gt.0.5d0*h.or.dabs(y).gt.0.5d0*h) then
 
 !rr=dsqrt(x*x+y*y)
 !cos=x/rr; sin=y/rr
! uu=rr*uuu*cos/radius; vv=rr*uuu*sin/radius 
!else
! uu=0.0d0; vv=0.0d0
!end if

!call tran3(rhoo,uu,vv,pp,ww, gg)
!c%value=ww
!c%gamma=gg
rhob=0.2228d0;  um=0.4d0; vm=0.0d0; pm=1.01325d0;  gammab=1.648d0
call tran3(rhob,um,vm,pm,wwb,gammab)
c%value=wwb
c%gamma=gammab

end function f_3


subroutine pick_pt(g_cell,pt)

implicit none
type(geo_info), intent(in) :: g_cell
real(8), dimension(2) :: pt

real(8) :: xy, alpha, beta

pt=error_data
select case(g_cell%g_type)
 case('xx ')
  pt(2)=0.0d0; xy=dfloat(g_cell%y_idx)
  xy=xy*h
  xy=dsqrt(radius*radius-xy*xy)
  if(g_cell%x_idx.lt.0) xy=-xy
  pt(1)=xy/h-dfloat(g_cell%x_idx)
 case('yy ')
  pt(1)=0.0d0; xy=dfloat(g_cell%x_idx)
  xy=xy*h
  xy=dsqrt(radius*radius-xy*xy)
  if(g_cell%y_idx.lt.0) xy=-xy
  pt(2)=xy/h-dfloat(g_cell%y_idx)
 case('xy ')
  alpha=dfloat(g_cell%x_idx); beta=dfloat(g_cell%y_idx)
  xy=dsqrt(alpha*alpha+beta*beta)
  alpha=alpha/xy; beta=beta/xy
  alpha=radius*alpha; beta=radius*beta
  pt(1)=alpha/h-dfloat(g_cell%x_idx)
  pt(2)=beta/h-dfloat(g_cell%y_idx)
 case default; call error_message
end select

end subroutine pick_pt


end module form_solution