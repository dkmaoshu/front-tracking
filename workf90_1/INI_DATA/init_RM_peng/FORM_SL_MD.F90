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
integer :: ii, jj, i, j, nn, first_idx, cases
real(8) :: pi
logical :: if_on_curve

real(8), parameter :: amplitude=0.1d0

public  set_in_smooth, set_discontinuity_curve

contains


subroutine set_in_smooth(cases)

implicit none

integer, intent(in) :: cases

real(8) :: gamma
real(8) :: rhot, ut, vt, pt, rhom, gamma_m, um, vm, pm,rhob, gamma_b
real(8), dimension(4) :: wwt, wwm, wwb

logical :: if_on_curve
real(8) :: ct, cb, cm, mu, m0, c_ast, dd
real(8), dimension(4) :: gt, gm
real(8) :: v_0, v_1, shock_speed, hc_check
integer :: i, j
type(state) :: ww

! real(8), dimension(2) :: nml
! character*12:: wtp
! type(state) :: l_nst, r_nst
! real(8) :: ve, rx,ux, vx, px

pi=4.0d0*datan(1.0d0)

select case(cases)

 case(1)

! See [Jacobs & Krivets], [Peng, Zabusky & Zhang] and [Courant & Friedrichs, pp149] for explanation of the following.

! The constant $\mu$.
  
!===========================================================
! For gamma=1.4
  !rhom=1.0d0; um=0.0d0; vm=0.0d0; pm=1.0d0; gamma_m=1.4d0    

 ! For gamma=1.2
 ! rhom=1.0d0; um=0.0d0; vm=0.0d0; pm=1.0d0; gamma_m=1.2d0

  ! ********Holmes(Air-SF6)********
!  rhom=0.95d0; um=0.0d0; vm=0.0d0; pm=0.8d0; gamma_m=1.4d0
  
  !******Holmes(Air-He)*********
  ! rhom=1.2d0; um=0.0d0; vm=0.0d0; pm=1.013d0; gamma_m=1.4d0

  ! ******Peng et al.(Air-SF6)******
  !rhom=1.18d0; um=0.0d0; vm=0.0d0; pm=1.0d0; gamma_m=1.4d0
  rhom=1.0d0; um=0.0d0; vm=0.0d0; pm=1.0d0; gamma_m=1.4d0


!===========================================================

  mu=dsqrt((gamma_m-1.0d0)/(gamma_m+1.0d0))

! The shock on the top is of Mach 5.
 ! m0=5.0 
 ! m0=1.2d0       !(Air-SF6)
 ! m0=1.52d0      !(Air-He)
  m0=1.3d0       !(Peng et al.)
  cm=dsqrt(gamma_m*pm/rhom)

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

  rhot=gamma_m*pt/(ct*ct)

  shock_speed=-v_0
  vt=v_1-v_0

  call tran3(rhom,um,vm,pm,wwm,gamma_m)  
  call tran3(rhot,ut,vt,pt,wwt,gamma_m)


! The following is the Hugoniot conditions check.

  do j=1, 4
   gt(j)=flux_in_y(j,wwt,gamma_m) 
   gm(j)=flux_in_y(j,wwm,gamma_m)
  end do
  hc_check=(gt(1)-gm(1))/(wwt(1)-wwm(1))
  hc_check=(gt(3)-gm(3))/(wwt(3)-wwm(3))
  hc_check=(gt(4)-gm(4))/(wwt(4)-wwm(4))

  wt%value=wwt; wm%value=wwm
!  wt%gamma=gamma_m; wm%gamma=gamma_m  

! The following sets the bottom side.

!  rhob=5.0d0 ; gamma_b=1.4d0 
!  rhob=4.84d0; gamma_b=1.09d0    !(Air-SF6)
!  rhob=0.167d0; gamma_b=1.63d0    !(Air-He)
!   rhob=5.2858d0; gamma_b=1.09d0    !(Peng et al.)
  rhob=4.52d0; gamma_b=1.4d0

!  rhob=4.52d0; gamma_b=1.4d0  !(Modified Peng et al.)
  call tran3(rhob,um,vm,pm,wwb,gamma_b)
  wb%value=wwb
 

! The following puts an upwards velocity to the setting.
!  vm=vm+3.12818921768967d0
 
!  vm=vm+0.21879607423073d0   !(Air-SF6)
  
!  vm=vm+1.13388459405120d0   !(Air-He)

  vm=vm+0.335973518044913d0    !(Peng et al.)
!  vm=vm+0.0d0    !(Peng et al.)

  call tran3(rhom,um,vm,pm,wwm,gamma_m)
  !vt=vt+3.12818921768967d0
  
!   vt=vt+0.21879607423073d0     !(Air-SF6)
 
!   vt=vt+1.13388459405120d0     !(Air-He)

  vt=vt+0.335973518044913d0      !(Peng et al.)
!  vt=vt+0.0d0      !(Peng et al.)

  call tran3(rhot,ut,vt,pt,wwt,gamma_m)
!  rhob=4.84d0   !(Density_SF6)
!  rhob=0.167d0  !(Density_He)
  rhob=4.52d0   !(Density_SF6, Peng et al.)
  call tran3(rhob,um,vm,pm,wwb,gamma_b)

  wt%value=wwt; wm%value=wwm; wb%value=wwb    
  wt%gamma=gamma_m; wm%gamma=gamma_m; wb%gamma=gamma_b
 
 
  
! The following is the initialization.
! First preset the smooth region that contains the shock.

!   wt=wm
   
!   nml=(/0.0d0,-1.0d0/)
!   call riemann(wt,wb,nml,3,0.0d0,0.0d0,l_nst,r_nst,ve,wtp)
!   call tran5(rx,ux,vx,px,l_nst%value)

  do i=0, nn
!   do j=-8*nn,8*nn
!   do j=-0.5*(30/4)*nn, 0.5*(30/4)*nn
!   do j=-9*nn, 9*nn/2
   do j=-nn-3, nn+3 ! For nonreflecting boundaries.
    dd=h*dfloat(j)
!    if(dd.gt.0.4d0+0.5d0*h) then
!    if(dd.gt.-5.0d0+0.5d0*h) then
!    if(dd.gt.6.0d0+0.5d0*h) then  !(Air/He)
	if(dd.gt.0.6d0+0.5d0*h) then  ! Peng et al.(Air/SF6)
     !uu(i,j)=wt
	 uu(i,j)%value=wt%value
	 uu(i,j)%gamma=gamma_m
    else
!	 if(dd.gt.0.4d0-0.5d0*h) then
!	 if(dd.lt.-5.0d0-0.5d0*h) then
!	 if(dd.lt.6.0d0-0.5d0*h) then    !(Air/He)
	 if(dd.lt.0.6d0-0.5d0*h) then  ! Peng et al.(Air/SF6)
	  !uu(i,j)=wm
	  uu(i,j)%value=wm%value
	  uu(i,j)%gamma=gamma_m           
     else
	  !uu(i,j)=0.5d0*(wm+wt)
	  uu(i,j)%value=0.5d0*(wwm+wwt)
	  uu(i,j)%gamma=gamma_m
     endif
    end if
   end do
  end do
  
! Then add the sinusoidal contact discontinuity to complete the intialization.     	 

  do i=0,nn
!   do j=-8*nn,8*nn
!    do j=-0.5*(30/4)*nn, 0.5*(30/4)*nn
!   do j=-9*nn, 9*nn/2
   do j=-nn-3, nn+3 ! For nonreflecting boundaries.
    call check_on_curve(i,j,if_on_curve)
    if(.not.if_on_curve) then
     dd=f_curve(dfloat(i)*h)
     if(dd.gt.dfloat(j)*h) then
      uu(i,j)%value=wwb; uu(i,j)%gamma=gamma_b
     end if
    end if
   end do
  end do

!  continue

 case(2)
  
 case default; print*, 'Something is wrong!!!'; pause
 
end select


contains


subroutine check_on_curve(x_idx,y_idx,if_on_curve)

implicit none
integer, intent(in) :: x_idx, y_idx
logical, intent(out) :: if_on_curve

real(8) :: x, y, rrp, rrm, dd

if_on_curve=.false.
x=dfloat(x_idx)*h; y=dfloat(y_idx)*h
rrp=f_curve(x+0.5d0*h)
rrm=f_curve(x-0.5d0*h)

if(rrm.lt.y+0.5d0*h.and.rrm.gt.y-0.5d0*h) then
 if_on_curve=.true.
 return
end if
if(rrp.lt.y+0.5d0*h.and.rrp.gt.y-0.5d0*h) then
 if_on_curve=.true.
 return
end if
 
end subroutine check_on_curve


end subroutine set_in_smooth


subroutine set_discontinuity_curve

implicit none

type(state) :: l_state, r_state, c1, c2
real(8), dimension(2) :: normal, pt1, pt2
real(8) :: velocity, posi_l, posi_r, posi_t, posi_b
real(8) :: xx, yy, tiny_x, tiny_y, hh, dd
real(8), dimension(4) :: sum
character*12 :: wave_type
character*8 :: cases
integer :: k, jj
type(state) :: fff
hh=h/dfloat(20)

call clean_up_grid_map(ggd_cell)
print*, 'Creat the discontinuity curve.'

cvv(1)%status='awake'
cvv(1)%cv_type='dblink'
cvv(1)%begin_end=0; cvv(1)%end_end=0
allocate(cvv(1)%begin)
nullifY(cvv(1)%begin%previous)
current=>cvv(1)%begin

k=1
! The first critical cell in the list is set be the leftest one.
do ii=nxll, nxll+nx-1
 posi_l=error_data; posi_r=error_data
 posi_t=error_data; posi_b=error_data
 call clean_up_g_cell(g_cell); call clean_up_g_cell(g_cell_1)
 call clean_up_p_cell(p_cell); call clean_up_p_cell(p_cell_1)

! Set geometric information.
 xx=dfloat(ii)*h; yy=f_curve(xx)
 if(yy.lt.0.0d0) then
  jj=idint(yy/h-0.5d0)
 else
  jj=idint(yy/h+0.5d0)
 end if
 yy=dfloat(jj)*h
 posi_l=(f_curve(xx-0.5d0*h)-yy)/h; posi_r=(f_curve(xx+0.5d0*h)-yy)/h

 cases='single  '

 if(posi_l.gt.0.5d0) then
  cases='two_up_l'; posi_l=posi_l-1.0d0
  if(xx.lt.1.0d0) then
   posi_t=(inv_f_curve_n(yy+0.5d0*h)-xx)/h
  else
   posi_t=(inv_f_curve_p(yy+0.5d0*h)-xx)/h 
  end if
 else
  if(posi_l.lt.-0.5d0) then
   cases='two_dn_l'; posi_l=posi_l+1.0d0
   if(xx.lt.1.0d0) then
    posi_b=(inv_f_curve_n(yy-0.5d0*h)-xx)/h
   else
    posi_b=(inv_f_curve_p(yy-0.5d0*h)-xx)/h 
   end if
  end if
 end if 
 
 if(posi_r.gt.0.5d0) then
  cases='two_up_r'; posi_r=posi_r-1.0d0
  if(xx.lt.1.0d0) then
   posi_t=(inv_f_curve_n(yy+0.5d0*h)-xx)/h
  else
   posi_t=(inv_f_curve_p(yy+0.5d0*h)-xx)/h 
  end if
 else
  if(posi_r.lt.-0.5d0) then
   cases='two_dn_r'; posi_r=posi_r+1.0d0
   if(xx.lt.1.0d0) then
    posi_b=(inv_f_curve_n(yy-0.5d0*h)-xx)/h
   else
    posi_b=(inv_f_curve_p(yy-0.5d0*h)-xx)/h 
   end if
  end if
 end if 
 
 select case(cases)

  case('single  ')
   g_cell%x_idx=ii; g_cell%y_idx=jj
   g_cell%edge(1)=4; g_cell%dis_posi(1)=posi_l
   g_cell%edge(2)=2; g_cell%dis_posi(2)=posi_r

  case('two_up_l')
   g_cell%x_idx=ii; g_cell%y_idx=jj+1
   g_cell%edge(1)=4; g_cell%dis_posi(1)=posi_l
   g_cell%edge(2)=1; g_cell%dis_posi(2)=posi_t
   g_cell_1%x_idx=ii; g_cell_1%y_idx=jj
   g_cell_1%edge(1)=3; g_cell_1%dis_posi(1)=posi_t
   g_cell_1%edge(2)=2; g_cell_1%dis_posi(2)=posi_r   

  case('two_dn_l')
   g_cell%x_idx=ii; g_cell%y_idx=jj-1
   g_cell%edge(1)=4; g_cell%dis_posi(1)=posi_l
   g_cell%edge(2)=3; g_cell%dis_posi(2)=posi_b
   g_cell_1%x_idx=ii; g_cell_1%y_idx=jj
   g_cell_1%edge(1)=1; g_cell_1%dis_posi(1)=posi_b
   g_cell_1%edge(2)=2; g_cell_1%dis_posi(2)=posi_r

  case('two_up_r')
   g_cell%x_idx=ii; g_cell%y_idx=jj
   g_cell%edge(1)=4; g_cell%dis_posi(1)=posi_l
   g_cell%edge(2)=3; g_cell%dis_posi(2)=posi_t
   g_cell_1%x_idx=ii; g_cell_1%y_idx=jj+1
   g_cell_1%edge(1)=1; g_cell_1%dis_posi(1)=posi_t
   g_cell_1%edge(2)=2; g_cell_1%dis_posi(2)=posi_r   

  case('two_dn_r')
   g_cell%x_idx=ii; g_cell%y_idx=jj
   g_cell%edge(1)=4; g_cell%dis_posi(1)=posi_l
   g_cell%edge(2)=1; g_cell%dis_posi(2)=posi_b
   g_cell_1%x_idx=ii; g_cell_1%y_idx=jj-1
   g_cell_1%edge(1)=3; g_cell_1%dis_posi(1)=posi_b
   g_cell_1%edge(2)=2; g_cell_1%dis_posi(2)=posi_r
   
  case default; call error_message     
         
 end select     

 call find_points_from_edges(g_cell)
 call find_type_from_edges(g_cell)
 if(cases.eq.'two_up_l'.or.cases.eq.'two_dn_l') then
  call find_points_from_edges(g_cell_1)
  call find_type_from_edges(g_cell_1)
 end if
 if(cases.eq.'two_up_r'.or.cases.eq.'two_dn_r') then
  call find_points_from_edges(g_cell_1)
  call find_type_from_edges(g_cell_1)
 end if

! Set physical information.

 fff=wm
 fff=wb
 p_cell%l_state=wm; p_cell%r_state=wb
 jj=g_cell%y_idx
 sum=0.0d0
 do i=-10, 10
  do j=-10, 10
   tiny_x=dfloat(i)*hh; tiny_y=dfloat(j)*hh
   xx=dfloat(ii)*h+tiny_x; yy=dfloat(jj)*h+tiny_y
   dd=f_curve(xx)
   if(dd.lt.yy) then
    sum=sum+wm%value
   else
    sum=sum+wb%value
   end if
  end do
 end do
 sum=sum/441.0d0
 
 p_cell%or_state%value=sum; p_cell%wv_nb=2
 p_cell%or_state%gamma=error_data


 if(cases.eq.'two_up_l'.or.cases.eq.'two_dn_l') then
  p_cell_1%l_state=wm; p_cell_1%r_state=wb
  jj=g_cell_1%y_idx; 

  sum=0.0d0
  do i=-10, 10
   do j=-10, 10
    tiny_x=dfloat(i)*hh; tiny_y=dfloat(j)*hh
    xx=dfloat(ii)*h+tiny_x; yy=dfloat(jj)*h+tiny_y
    dd=f_curve(xx)
    if(dd.lt.yy) then
     sum=sum+wm%value
    else
     sum=sum+wb%value
    end if
   end do
  end do
  sum=sum/441.0d0

  p_cell_1%or_state%value=sum; p_cell_1%wv_nb=2
  p_cell_1%or_state%gamma=error_data
 
 end if

 if(cases.eq.'two_up_r'.or.cases.eq.'two_dn_r') then
  p_cell_1%l_state=wm; p_cell_1%r_state=wb
  jj=g_cell_1%y_idx
  
  sum=0.0d0
  do i=-10, 10
   do j=-10, 10
    tiny_x=dfloat(i)*hh; tiny_y=dfloat(j)*hh
    xx=dfloat(ii)*h+tiny_x; yy=dfloat(jj)*h+tiny_y
    dd=f_curve(xx)
    if(dd.lt.yy) then
     sum=sum+wm%value
    else
     sum=sum+wb%value
    end if
   end do
  end do
  sum=sum/441.0d0

  p_cell_1%or_state%value=sum; p_cell_1%wv_nb=2 
  p_cell_1%or_state%gamma=error_data 
  
 end if

 allocate(temp); nullify(temp%previous); nullify(temp%next)
 temp%g_cell=g_cell; temp%p_cell=p_cell
 temp%address%idx=k; k=k+1
 temp%address%cv_nmb=1
 temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)

 jj=g_cell%y_idx
 uu(ii,jj)=error_data
 ggd_cell(ii,jj)%region='crit'
 do j=1,4
  ggd_cell(ii,jj)%ccpt(j)%address=temp%address
  ggd_cell(ii,jj)%ccpt(j)%side=g_cell%point(j)
 end do

 temp%previous=>current
 current%next=>temp
 current=>temp

 if(cases.eq.'two_up_l'.or.cases.eq.'two_dn_l') then
  allocate(temp); nullify(temp%previous); nullify(temp%next)
  temp%g_cell=g_cell_1; temp%p_cell=p_cell_1
  temp%address%idx=k; k=k+1
  temp%address%cv_nmb=1
  temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)
  
  jj=g_cell_1%y_idx
  uu(ii,jj)=error_data
  ggd_cell(ii,jj)%region='crit'
  do j=1,4
   ggd_cell(ii,jj)%ccpt(j)%address=temp%address
   ggd_cell(ii,jj)%ccpt(j)%side=g_cell_1%point(j)
  end do

  temp%previous=>current
  current%next=>temp
  current=>temp
 end if

 if(cases.eq.'two_up_r'.or.cases.eq.'two_dn_r') then
  allocate(temp); nullify(temP%previous); nullify(temp%next)
  temp%g_cell=g_cell_1; temp%p_cell=p_cell_1
  temp%address%idx=k; k=k+1
  temp%address%cv_nmb=1
  temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)
  
  jj=g_cell_1%y_idx
  uu(ii,jj)=error_data
  ggd_cell(ii,jj)%region='crit'
  do j=1,4
   ggd_cell(ii,jj)%ccpt(j)%address=temp%address
   ggd_cell(ii,jj)%ccpt(j)%side=g_cell_1%point(j)
  end do

  temp%previous=>current
  current%next=>temp
  current=>temp
 end if

end do

nullify(cvv(1)%begin%next%previous)
allocate(cvv(1)%eend)	
nullify(cvv(1)%eend%next)
cvv(1)%eend%previous=>current
cvv(1)%total=k-1
cvv(1)%wave=2

print*, 'The discontinuity curve has been created.' 
! Created the 'curve_number'th discontinuity curve.

end subroutine set_discontinuity_curve


function f_curve(x) result(y)

implicit none

real(8), intent(in) :: x
real(8) :: y

!y=0.1d0*dcos(2.0d0*pi*(x-0.5d0))

    !***** Air/SF6  !**  [y=a(t)cos((2*pi/lamda)*(x-lamda/2))+(-(y_width/2)+interface_from_x-axis)]
!y=0.24d0*dcos(((2.0d0*pi)/3.75d0)*(x-1.875d0))-5.75d0
    !****** Air/He  !******
!y=0.2d0*dcos(((2.0d0*pi)/4.0d0)*(x-2.0d0))+4.9d0

!  !****Peng et al* Air/SF6  !**
!y=0.1d0*dcos(((2.0d0*pi)/2.0d0)*(x-1.0d0))

  !****Second Peng et al* Air/SF6  !**
y=amplitude*dcos(((2.0d0*pi)/2.0d0)*(x-1.0d0))

end function f_curve


function inv_f_curve_n(y) result(x)

implicit none

real(8), intent(in) :: y
real(8) :: x

!x=-dacos(y/0.1d0)/pi/2.0d0+0.5d0
!x=-dacos(y/0.24d0)/1.675d0+0.5d0
  !***** Air-SF6  !******
!x=-3.75d0*dacos((y+5.75d0)/0.24d0)/pi/2.0d0+1.875d0
   !***** Air-He  !****
! x=-4.0d0*dacos((y-4.90d0)/0.2d0)/pi/2.0d0+2.0d0

!!****Peng et al* Air/SF6  !**
!x=-2.0d0*dacos(y/0.1d0)/pi/2.0d0+1.0d0

!****Second Peng et al* Air/SF6  !**
x=-2.0d0*dacos(y/amplitude)/pi/2.0d0+1.0d0

end function inv_f_curve_n


function inv_f_curve_p(y) result(x)

implicit none

real(8), intent(in) :: y
real(8) :: x

!x=dacos(y/0.1d0)/pi/2.0d0+0.5d0
!x=dacos((y+5.75d0)/0.24d0)*1.0d0/1.675d0+1.675d0
  !***** Air-SF6  !*****
!x=3.75d0*dacos((y+5.75d0)/0.24d0)/pi/2.0d0+1.875d0
  !***** Air-He !*****
! x=4.0d0*dacos((y-4.9d0)/0.2d0)/pi/2.0d0+2.0d0

!!****Peng et al* Air/SF6  !**
!x=2.0d0*dacos(y/0.1d0)/pi/2.0d0+1.0d0

!****Second Peng et al* Air/SF6  !**
x=2.0d0*dacos(y/amplitude)/pi/2.0d0+1.0d0


end function inv_f_curve_p


end module form_solution