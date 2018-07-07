module output_show_HR

use solution
! 'solution.f90'

implicit none

integer, parameter :: res_numb=2

public show_out_hr


contains


subroutine show_out_hr
! This subroutine outputs the numerical solution with high resolution for display by
! Matlab.

implicit none

type(state), dimension(:,:), allocatable :: uu_show
real(8), dimension(:,:), allocatable :: rho_show, u_show, v_show, p_show
type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
!type(node_info) :: n_cell

integer:: head_mark, i, j, ix, jy, ii, jj, single
logical :: head_switch
integer :: left_bound, right_bound, bottom, top, local_center_x, local_center_y
real(8), dimension(2) :: local_point
type(state), dimension(-res_numb:res_numb,-res_numb:res_numb) ::  &
             local_u, local_u_left, local_u_right
real(8) :: length
type(state) :: x_slopl,x_slopr,y_slopu,y_slopd
character*8 :: side

type(state), dimension(-3:3) :: ux, uy
real(8), dimension(4) :: value
real(8) :: hh=1.0d0/(2.0d0*dfloat(res_numb)+1.0d0)

! The numerical solution is output on a grid 10 times finer than the original grid. Each 
! grid cell of the original grid is fined with 10 grid cells in both $x$- and $y$-direction.
! Thus, in the numerical solution in the critical cell cells are ouput according to whether
! the finer grid cells are either on the left or right side of the tracked front.

! First, determine the size of the finer grid.
left_bound=(2*res_numb+1)*nxll-res_numb 
right_bound=(2*res_numb+1)*(nxll+nx-1)+res_numb
bottom=(2*res_numb+1)*nyll-res_numb
top=(2*res_numb+1)*(nyll+ny-1)+res_numb
allocate(uu_show(left_bound:right_bound,bottom:top))
allocate(rho_show(left_bound:right_bound,bottom:top))
allocate(u_show(left_bound:right_bound,bottom:top))
allocate(v_show(left_bound:right_bound,bottom:top))
allocate(p_show(left_bound:right_bound,bottom:top))

do ii=left_bound,right_bound
 do jj=bottom,top
  uu_show(ii,jj)=error_data
 end do
end do

do i=nxll,nxll+nx-1
  do j=nyll,nyll+ny-1
  if(ggd_cell(i,j)%region.eq.'smth') then

! The following  computes the left and right pieces of the smooth solution.
   call smth_dpr(i,j,temp,ux,'xx ','      ','ordinary',.true.,3)
   call smth_dpr(i,j,temp,uy,'yy ','      ','ordinary',.true.,3)
   x_slopl=ux(0)-ux(-1); y_slopd=uy(0)-uy(-1)
   x_slopr=ux(1)-ux(0); y_slopu=uy(1)-uy(0)
   do ii=-res_numb,res_numb
    do jj=-res_numb,res_numb
     local_point(1)=dfloat(ii)*hh; local_point(2)=dfloat(jj)*hh
     if(ii.lt.0) then
	  if(jj.lt.0) then
       local_u(ii,jj)=uu(i,j)+local_point(1)*x_slopl+local_point(2)*y_slopd
      else
       local_u(ii,jj)=uu(i,j)+local_point(1)*x_slopl+local_point(2)*y_slopu
      end if
     else
	  if(jj.lt.0) then
       local_u(ii,jj)=uu(i,j)+local_point(1)*x_slopr+local_point(2)*y_slopd
      else
       local_u(ii,jj)=uu(i,j)+local_point(1)*x_slopr+local_point(2)*y_slopu
      end if
     end if
    end do
   end do
   local_center_x=i*(2*res_numb+1); local_center_y=(2*res_numb+1)*j
   uu_show(local_center_x-res_numb:local_center_x+res_numb, &
           local_center_y-res_numb:local_center_y+res_numb)=local_u
  end if 
 end do
end do

! call check_list_c(2,'down  ')
! call check_ring('counter ')

do i=1,curves_number
 if(cvv(i)%status.eq.'awake') then
  nullify(temp)
   
!   sum=0.d0

!   open(11,file='d:\workf90_1\xxx.dat')

  if(associated(cvv(i)%begin%next)) then
   temp=>cvv(i)%begin%next
   head_mark=temp%address%idx; head_switch=.true.
  end if
  do while(associated(temp).and.head_switch) 
   g_cell=temp%g_cell; p_cell=temp%p_cell

!	write(11,'(2i5,f16.9)') g_cell%x_idx, g_cell%y_idx, p_cell%or_state
!	write(11,'(f16.9)') p_cell%or_state%value
!	write(11,*)
    
   ix=g_cell%x_idx; jy=g_cell%y_idx

! First, compute the left and right pieces of the smooth solution.
   call smth_dpr(i,j,temp,ux,'xx ','left  ','critical',.true.,3)
   call smth_dpr(i,j,temp,uy,'yy ','left  ','critical',.true.,3)
   x_slopl=ux(0)-ux(-1); y_slopd=uy(0)-uy(-1)
   x_slopr=ux(1)-ux(0); y_slopu=uy(1)-uy(0)
   do ii=-res_numb,res_numb
    do jj=-res_numb,res_numb
     local_point(1)=dfloat(ii)*hh; local_point(2)=dfloat(jj)*hh
     if(ii.lt.0) then
	  if(jj.lt.0) then
       local_u_left(ii,jj)=ux(0)+local_point(1)*x_slopl+local_point(2)*y_slopd
      else
       local_u_left(ii,jj)=ux(0)+local_point(1)*x_slopl+local_point(2)*y_slopu
      end if
     else
	  if(jj.lt.0) then
       local_u_left(ii,jj)=ux(0)+local_point(1)*x_slopr+local_point(2)*y_slopd
      else
       local_u_left(ii,jj)=ux(0)+local_point(1)*x_slopr+local_point(2)*y_slopu
      end if
     end if
    end do
   end do
   call smth_dpr(i,j,temp,ux,'xx ','right ','critical',.true.,3)
   call smth_dpr(i,j,temp,uy,'yy ','right ','critical',.true.,3)
   x_slopl=ux(0)-ux(-1); y_slopd=uy(0)-uy(-1)
   x_slopr=ux(1)-ux(0); y_slopu=uy(1)-uy(0)
   do ii=-res_numb,res_numb
    do jj=-res_numb,res_numb
     local_point(1)=dfloat(ii)*hh; local_point(2)=dfloat(jj)*hh
     if(ii.lt.0) then
	  if(jj.lt.0) then
       local_u_right(ii,jj)=ux(0)+local_point(1)*x_slopl+local_point(2)*y_slopd
      else
       local_u_right(ii,jj)=ux(0)+local_point(1)*x_slopl+local_point(2)*y_slopu
      end if
     else
	  if(jj.lt.0) then
       local_u_right(ii,jj)=ux(0)+local_point(1)*x_slopr+local_point(2)*y_slopd
      else
       local_u_right(ii,jj)=ux(0)+local_point(1)*x_slopr+local_point(2)*y_slopu
      end if
     end if
    end do
   end do

!   call compute_curve_length(g_cell,length)

!   if(length.gt.0.02d0) then
   local_center_x=(2*res_numb+1)*ix; local_center_y=(2*res_numb+1)*jy
   do ii=-res_numb,res_numb
    do jj=-res_numb,res_numb
     local_point(1)=dfloat(ii)*hh; local_point(2)=dfloat(jj)*hh
	 call check_region(local_point,side)
!     call check_left_right(local_point,side)
     select case(side)
      case('left    '); local_u(ii,jj)=local_u_left(ii,jj)
      case('right   '); local_u(ii,jj)=local_u_right(ii,jj)
	  case('neither '); local_u(ii,jj)=error_data
      case default; call error_message
     end select
    end do
   end do
 !  else
 !   if(g_cell%g_type.ne.'xy ') call error_message
 !   call find_single(g_cell,single)
 !   select case(g_cell%point(single))
 !    case('left  '); local_u=local_u_right
 !    case('right '); local_u=local_u_left
 !    case default; call error_message
 !   end select
 !  end if
   do ii=-res_numb,res_numb
    do jj=-res_numb,res_numb
     if(local_u(ii,jj).gt.0.9d0*error_data) then
	  uu_show(local_center_x+ii,local_center_y+jj)=local_u(ii,jj)
     end if
    end do
   end do   
!   local_center_x=(2*res_numb+1)*ix; local_center_y=(2*res_numb+1)*jy
!   uu_show(local_center_x-res_numb:local_center_x+res_numb,&
!           local_center_y-res_numb:local_center_y+res_numb)=local_u
    
   temp=>temp%next
   if(associated(temp)) then
    head_switch=(temp%address%idx.ne.head_mark)
   else
    exit
   end if
  end do
  
!   close(11)

 end if
end do

do i=left_bound, right_bound
 do j=bottom,top
  value=uu_show(i,j)%value
  rho_show(i,j)=value(1)
  u_show(i,j)=uf(value)
  v_show(i,j)=vf(value)
  p_show(i,j)=pf(value)
 end do
end do

! call scanning(ggd_cell)

!do i=1,nodes_number
! if(ndd(i)%status.eq.'awake') then
!  n_cell=ndd(i)%n_cell
!  u_show(n_cell%x_idx,n_cell%y_idx)=n_cell%or_state%value
! end if
!end do

open(2, file='d:\workf90_1\show\showrho_hr.dat')
do ii=left_bound, right_bound
 do jj=bottom, top
  write(2,'(f12.6)') rho_show(ii,jj)
 end do
end do  
close(2)

open(2, file='d:\workf90_1\show\showu_hr.dat')
do ii=left_bound, right_bound
 do jj=bottom, top
  write(2,'(f12.6)') u_show(ii,jj)
 end do
end do  
close(2)

open(2, file='d:\workf90_1\show\showv_hr.dat')
do ii=left_bound, right_bound
 do jj=bottom, top
  write(2,'(f12.6)') v_show(ii,jj)
 end do
end do  
close(2)

open(2, file='d:\workf90_1\show\showp_hr.dat')
do ii=left_bound, right_bound
 do jj=bottom, top
  write(2,'(f12.6)') p_show(ii,jj)
 end do
end do  
close(2)

open(2, file='d:\workf90_1\show\showres_hr.dat')
write(2,'(i5)') res_numb
close(2)

! Ouput for show of solution.

!call check_list(cvv(1),'down  ')

!deallocate(uu); deallocate(cvv)
!deallocate(ggd_cell); deallocate(ndd)


contains


subroutine check_region(local_point,side)

implicit none

real(8), dimension(2) :: local_point
character*8, intent(out) :: side

type(critical_cell), pointer :: temp_g
type(geo_info) :: g_cell_g
type(adss_info) :: adss_inf
character*8 :: side_1, side_2

call check_left_right(local_point,g_cell,side)

select case(side)

 case('left    ')
  adss_inf=temp%l_stk
  if(adss_inf%cv_nmb.ne.0) then
   call visit(adss_inf,temp_g)
   g_cell_g=temp_g%g_cell
  else
   return
  end if

 case('right   ')
  adss_inf=temp%r_stk
  if(adss_inf%cv_nmb.ne.0) then
   call visit(adss_inf,temp_g)
   g_cell_g=temp_g%g_cell
  else
   return
  end if   

 case default; call error_message
 
end select
   
call check_left_right(local_point,g_cell_g,side_1)
if(temp_g%l_stk.eq.temp%address) then
 side_2='left    ' 
else
 if(temp_g%r_stk.eq.temp%address) then
  side_2='right   '
 else
  call error_message
 end if
end if

if(side_1.ne.side_2) side='neither'

end subroutine check_region


subroutine check_left_right(local_point,g_cell,side)

implicit none

real(8), dimension(2), intent(in) :: local_point
type(geo_info), intent(in) :: g_cell
character*8, intent(out) :: side

real(8) :: aa

real(8), dimension(2) :: pt1, pt2, vec_1, vec_2

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)

vec_1(1)=pt2(1)-pt1(1); vec_1(2)=pt2(2)-pt1(2)
vec_2(1)=local_point(1)-pt1(1); vec_2(2)=local_point(2)-pt1(2)
aa=outer_product(vec_1,vec_2)
if(aa.gt.0.0d0) then
 side='left    '
else
 side='right   '
end if

end subroutine check_left_right


function outer_product(vec_1,vec_2) result(c)

implicit none
real(8), dimension(2), intent(in) :: vec_1, vec_2
real(8) :: c

c=vec_1(1)*vec_2(2)-vec_2(1)*vec_1(2)

end function outer_product


end subroutine show_out_hr


end module output_show_HR