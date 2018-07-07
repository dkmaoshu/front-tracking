module produce_polynomial_4_xx_yy
! This module is for producing polynomial $p(x,y)$ on the left and right sides of a critical
! cell.

use interpolation
! 'interpl.f90'

use solu_comput
! 'solu_com.f90'

implicit none


public  p_side, lagrange_state
private lagrange, lagrange_derivative


contains


function p_side(temp,x,y,order,side) result(c)
! The polynomial function associated with the critical cell pointed by 'temp'. Symbolics in the
! remarks here can be refered to 'reset.tex'. 

implicit none

type(auxi_crit_cell), pointer :: temp
! The pointer pointing to the critical cell concerned.
real(8), intent(in) :: x, y
! The two variables.
integer, intent(in) :: order
! The order of the polynomial.
character*6, intent(in) :: side
! The side on which the polynomial is constructed.
type(state) :: c

type(auxi_crit_cell), pointer :: temp_nb
integer, dimension(-1:1) :: i0, j0
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(state), dimension(-3:3) :: uxy
type(state), dimension(-1:1) :: uuu
character*3 :: xy_dir
integer :: i, ii
type(state), dimension(-1:1) :: uuxy
real(8), dimension(-1:1) :: xxyy
real(8) :: xy
type(state) :: xxx

g_cell=temp%g_cell; p_cell=temp%p_cell
if(g_cell%g_type.ne.'xx'.and.g_cell%g_type.ne.'yy') call error_message
! The function is valid only for 'xx'- or 'yy'-type cirtical cells.

! Preparation.
xxyy(-1)=-1.0d0; xxyy(0)=0.0d0; xxyy(1)=1.0d0
xy_dir=g_cell%g_type
do ii=-1,1
 uuxy(ii)%value=error_data
end do
select case(xy_dir)
 case('xx '); xy=x
 case('yy '); xy=y
end select

! Produce reconstruction polynomial $p^{n,\pm}_j(x)$.
i0(0)=g_cell%x_idx; j0(0)=g_cell%y_idx
call smth_dpr(i0(0),j0(0),temp,uxy,xy_dir,side,'critical',.true.,order,'full    ')
do ii=-1,1
 uuu(ii)=uxy(ii)
end do
uuxy(0)=state_inter(uuu,xxyy,xy)

! Produce reconstruction polynomials $p^{n,\pm}_{j+1}(x)$ and $p^{n,\pm}_{j-1}(x)$.
temp_nb=>temp%p_nxt_dr
if(associated(temp_nb)) then
 g_cell=temp_nb%g_cell; p_cell=temp_nb%p_cell
 i0(-1)=g_cell%x_idx; j0(-1)=g_cell%y_idx
 call smth_dpr(i0(-1),j0(-1),temp_nb,uxy,xy_dir,side,'critical',.true.,order,'full    ')
 select case(xy_dir)
  case('xx ') 
   if(i0(-1).gt.i0(0)) then
    xy=xy-1.0d0
   else
    if(i0(-1).lt.i0(0)) then
     xy=xy+1.0d0
    end if
   end if
  case('yy ')
   if(j0(-1).gt.j0(0)) then
    xy=xy-1.0d0
   else
    if(j0(-1).lt.j0(0)) then
     xy=xy+1.0d0
    end if
   end if
 end select
 do ii=-1,1
  uuu(ii)=uxy(ii)
 end do
 uuxy(-1)=state_inter(uuu,xxyy,xy)
end if

select case(xy_dir)
 case('xx '); xy=x
 case('yy '); xy=y
end select
temp_nb=>temp%n_nxt_dr
if(associated(temp_nb)) then
 g_cell=temp_nb%g_cell; p_cell=temp_nb%p_cell
 i0(1)=g_cell%x_idx; j0(1)=g_cell%y_idx
 call smth_dpr(i0(1),j0(1),temp_nb,uxy,xy_dir,side,'critical',.true.,order,'full    ')
 select case(xy_dir)
  case('xx ') 
   if(i0(1).gt.i0(0)) then
    xy=xy-1.0d0
   else
    if(i0(1).lt.i0(0)) then
     xy=xy+1.0d0
    end if
   end if
  case('yy ')
   if(j0(1).gt.j0(0)) then
    xy=xy-1.0d0
   else
    if(j0(1).lt.j0(0)) then
     xy=xy+1.0d0
    end if
   end if
 end select
 do ii=-1,1 
  uuu(ii)=uxy(ii)
 end do
 uuxy(1)=state_inter(uuu,xxyy,xy)
end if

! Produce $p^{n,-}(x,y)
select case(xy_dir)
 case('yy ')
  do i=-1,1 
   xxyy(i)=dfloat(i0(i)-i0(0))
  end do
  xy=x
 case('xx ')
  do i=-1,1 
   xxyy(i)=dfloat(j0(i)-j0(0))
  end do
  xy=y
end select

xxx=state_inter(uuxy,xxyy,xy)
  
c=xxx


contains


function state_inter(uuu,xxyy,xy) result(c)
! The interpolation function via prime function on physical states.

implicit none

type(state), dimension(-1:1), intent(in) :: uuu
real(8), dimension(-1:1), intent(in) :: xxyy
real(8), intent(in) :: xy
type(state) :: c

integer :: i, j
real(8), dimension(0:2) :: vv, xyp
real(8), dimension(0:1) :: vvv, xypp
character*8 :: direction
real(8) :: cc
character*14 :: cases

if(xxyy(0).gt.xxyy(-1)) then
 direction='upward  '
else
 direction='downward'
end if

cases='no_missing    '
if(uuu(-1).lt.0.9d0*error_data) then
 if(uuu(1).lt.0.9d0*error_data) then
  cases='both_missing  '
 else
  cases='left_missing  '
 end if
else
 if(uuu(1).lt.0.9d0*error_data) then
  cases='right_missing '
 end if
end if        

select case(cases)
 
 case('no_missing    ')
  do j=1,4
   do i=-1,1
    vv(i+1)=uuu(i)%value(j); xyp(i+1)=xxyy(i)
   end do
   cc=interpol_via_prime_fun(vv,xyp,xy,2,0.5d0,direction)
   c%value(j)=cc
   c%gamma=uuu(0)%gamma
  end do

 case('left_missing  ')
  do j=1,4
   do i=0,1
    vvv(i)=uuu(i)%value(j); xypp(i)=xxyy(i)
   end do
   cc=interpol_via_prime_fun(vvv,xypp,xy,1,0.5d0,direction)
   c%value(j)=cc
   c%gamma=uuu(0)%gamma
  end do

 case('right_missing ')
  do j=1,4
   do i=-1,0
    vvv(i+1)=uuu(i)%value(j); xypp(i+1)=xxyy(i)
   end do
   cc=interpol_via_prime_fun(vvv,xypp,xy,1,0.5d0,direction)
   c%value(j)=cc
   c%gamma=uuu(0)%gamma
  end do
  
 case('both_missing  ')
  c=uuu(0)

end select      

end function state_inter


end function p_side


function lagrange_state(u,x,y,n) result(c)
! Lagrange interpolation applied to each component of physical states.

implicit none
integer, intent(in) :: n
type(state), dimension(0:n), intent(in) :: u
real(8), dimension(0:n), intent(in) :: x
real(8), intent(in) :: y
type(state) :: c

integer :: i, j
real(8), dimension(0:n) :: vv
real(8) :: cc

do j=1,4
 do i=0,n
  vv(i)=u(i)%value(j)
 end do
 cc=lagrange(vv,x,y,n)
 c%value(j)=cc
 c%gamma=u(0)%gamma
end do

end function lagrange_state


end module produce_polynomial_4_xx_yy



