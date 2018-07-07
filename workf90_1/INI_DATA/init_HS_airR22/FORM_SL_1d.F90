module form_solution_1d
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

public  set_in_smooth, set_discontinuity_curve, cases


contains


subroutine set_in_smooth

implicit none

real(8) :: rhot, ut, vt, pt, rhom, um, vm, pm, rhob, gammam, gammab, gammat
real(8), dimension(4) :: wwt, wwm, wwb

real(8) :: ct, cb, cm, mu, m0, c_ast, dd
real(8), dimension(4) :: gt, gm
real(8) :: v_0, v_1, shock_speed, hc_check, xx0, xx1
integer :: i, j
type(state) :: ww


rhom=1.225d0;  um=0.0d0;       vm=0.0d0; pm=1.01325d0; gammam=1.4d0
! rhot=1.6861d0; ut=-0.494137d0; vt=0.0d0; pt=2.50638d0; gammat=1.4d0
rhob=3.8635d0;                                         gammab=1.249d0
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
  xx0=(dfloat(i)-0.5d0)*h; xx1=(dfloat(i)+0.5d0)*h
  if(xx0.gt.-radius.and.xx1.lt.radius) then
   uu(i,j)%value=wb%value
   uu(i,j)%gamma=gammab
  else
   if(xx1.lt.-radius.or.xx0.gt.radius) then
    if(xx0.gt.0.3d0) then
     uu(i,j)%value=wwt
     uu(i,j)%gamma=gammat
    else
     if(xx0.lt.0.3d0.and.xx1.gt.0.3d0) then  
      uu(i,j)%value=0.5d0*(wwm+wwt)
      uu(i,j)%gamma=gammam
     else 
      uu(i,j)%value=wwm
      uu(i,j)%gamma=gammam 
     end if
    end if
   end if
  end if 
 end do
end do

end subroutine set_in_smooth


subroutine set_discontinuity_curve

implicit none

type(state) :: l_state, r_state, c1, c2
real(8), dimension(2) :: normal, pt1, pt2
real(8) :: velocity
character*12 :: wave_type
integer :: ii0, j

real(8) :: xposi, yposi, zzz

call clean_up_grid_map(ggd_cell)

zzz=radius/h
ii0=int(radius/h+0.5d0)
zzz=zzz-dfloat(ii0)

print*, 'Creat the first discontinuity curve.'

cvv(1)%status='awake'
cvv(1)%cv_type='dblink  '
cvv(1)%begin_end=0; cvv(1)%end_end=0
allocate(cvv(1)%begin)
nullifY(cvv(1)%begin%previous)
current=>cvv(1)%begin

do jj=-5, 5

! Geometry information.
 g_cell%x_idx=-ii0; g_cell%y_idx=jj
 yy=-0.5d0*h
 g_cell%edge(1)=1; g_cell%edge(2)=3
 g_cell%dis_posi(1)=1.0d0+zzz; g_cell%dis_posi(2)=1.0d0+zzz
 call find_points_from_edges(g_cell)
 call find_type_from_edges(g_cell)

! Physical information.
 p_cell%l_state=wm; p_cell%r_state=wb
 p_cell%or_state=(zzz+0.5d0)*wm+(0.5d0-zzz)*wb; p_cell%or_state%gamma=error_data
 p_cell%wv_nb=2

! Link and map information.
 allocate(temp)
 temp%g_cell=g_cell; temp%p_cell=p_cell
 temp%address%idx=jj+6
 temp%address%cv_nmb=1
 temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)
 uu(-ii0,jj)=error_data
 ggd_cell(-ii0,jj)%region='crit'
 do j=1,4
  ggd_cell(-ii0,jj)%ccpt(j)%address=temp%address
  ggd_cell(-ii0,jj)%ccpt(j)%side=g_cell%point(j)
 end do

! Linking
 temp%previous=>current
 current%next=>temp
 current=>temp

end do

nullify(cvv(1)%begin%next%previous)
allocate(cvv(1)%eend)	
nullify(cvv(1)%eend%next)
!cvv(1)%eend%previous=>current
cvv(1)%total=11
cvv(1)%wave=2

print*, 'The first discontinuity curve has been created.' 

print*, ' '

print*, 'Creat the second discontinuity curve.'

cvv(2)%status='awake'
cvv(2)%cv_type='dblink  '
cvv(2)%begin_end=0; cvv(1)%end_end=0
allocate(cvv(2)%begin)
nullifY(cvv(2)%begin%previous)
current=>cvv(2)%begin

do jj=-5, 5

! Geometry information.
 g_cell%x_idx=ii0; g_cell%y_idx=jj
 yy=-0.5d0*h
 g_cell%edge(1)=1; g_cell%edge(2)=3
 g_cell%dis_posi(1)=zzz; g_cell%dis_posi(2)=zzz
 call find_points_from_edges(g_cell)
 call find_type_from_edges(g_cell)

! Physical information.
 p_cell%l_state=wb; p_cell%r_state=wm
 p_cell%or_state=(zzz+0.5d0)*wm+(0.5d0-zzz)*wb; p_cell%or_state%gamma=error_data
 p_cell%wv_nb=2

! Link and map information.
 allocate(temp)
 temp%g_cell=g_cell; temp%p_cell=p_cell
 temp%address%idx=jj+6
 temp%address%cv_nmb=2
 temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)
 uu(ii0,jj)=error_data
 ggd_cell(ii0,jj)%region='crit'
 do j=1,4
  ggd_cell(ii0,jj)%ccpt(j)%address=temp%address
  ggd_cell(ii0,jj)%ccpt(j)%side=g_cell%point(j)
 end do

! Linking
 temp%previous=>current
 current%next=>temp
 current=>temp

end do

nullify(cvv(2)%begin%next%previous)
allocate(cvv(2)%eend)	
nullify(cvv(2)%eend%next)
!cvv(1)%eend%previous=>current
cvv(2)%total=11
cvv(2)%wave=2

print*, 'The second discontinuity curve has been created.' 

end subroutine set_discontinuity_curve


end module form_solution_1d