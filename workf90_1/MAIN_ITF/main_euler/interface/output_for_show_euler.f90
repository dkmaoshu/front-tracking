module output_show

use solution
! 'solution.f90'

use physical_state
! 'eulernew.f90'

use scan_conservation_on_list
! 'sc_conls.f90'
implicit none

real(8), dimension(:,:), allocatable :: rho_show, u_show, v_show, p_show, vorticity_show
real(8), dimension(:,:), allocatable :: mach_number, acoustic_impedance

public  show_out
private scan_shw, curve_files, compute_sound_speed


contains


subroutine show_out
! This subroutine outputs the numerical solution for displayed by
! Matlab.

implicit none

type(state), dimension(:,:), allocatable :: uu_show
!type(state), dimension(-120:120,-120:120) :: uu_show

real(8), dimension(:,:), allocatable :: velo_radial, velo_tangential

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(node_info) :: n_cell

integer:: head_mark, i, j, ix, jy, l, lc, ll, m
real(8), dimension(0:10000) :: x, y
character*38 :: file(10)
logical :: head_switch
real(8), dimension(4) :: value

 real(8) :: xx, yy, rr, aaa, bbb
 real(8) :: u_radial_l, u_radial_r, u_radial, u_tangen_l, u_tangen_r, u_tangen
 type(state) :: xxx, yyy
 integer :: i0, j0

! real(8) :: xxxx, yyyy

! real(8) :: sum

! Output for the of numerical solution for show in global.

allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(cvn))
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ndd(ndn))
call give_sth(0)
call give_cv
call give_gcl
call give_nd

allocate(uu_show(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(rho_show(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(u_show(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(v_show(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(p_show(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(velo_radial(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(velo_tangential(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(vorticity_show(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(mach_number(nxll:nxll+nx-1,nyll:nyll+ny-1))
allocate(acoustic_impedance(nxll:nxll+nx-1,nyll:nyll+ny-1))
  
!call hidden_state_check_c
   
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  uu_show(i,j)=error_data
 end do
end do
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  if(ggd_cell(i,j)%region.eq.'smth') uu_show(i,j)=uu(i,j)
 end do
end do

! call check_list_c(1,'down  ')
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
    
!   write(11,'(2i5,f16.9)') g_cell%x_idx, g_cell%y_idx, p_cell%or_state9
!   write(11,'(f16.9)') p_cell%or_state%value
!   write(11,*)
     	 
    
   ix=g_cell%x_idx; jy=g_cell%y_idx
   if(uu_show(ix,jy).lt.0.9d0*error_data) then
    uu_show(ix,jy)=p_cell%or_state
     
!     uu_show(ix,jy)=p_cell%r_state	 
!     uu_show(ix,jy)=0.5d0*(uu_show(ix,jy)+p_cell%l_state)
!     sum=sum+u_show(ix,jy)*h
    
   else
    uu_show(ix,jy)=uu_show(ix,jy)+p_cell%or_state
    if(temp%l_stk%cv_nmb.le.0.and.temp%r_stk%cv_nmb.le.0) then
     call error_message
    end if 
   end if
   if(temp%l_stk%cv_nmb.gt.0) then
    uu_show(ix,jy)=uu_show(ix,jy)-0.5d0*p_cell%l_state        	 
   end if
   if(temp%r_stk%cv_nmb.gt.0) then
    uu_show(ix,jy)=uu_show(ix,jy)-0.5d0*p_cell%r_state        	 
   end if
    
    xx=dfloat(ix)*h; yy=dfloat(jy)*h
    rr=dsqrt(xx*xx+yy*yy)
	xx=xx/rr; yy=yy/rr
    u_radial_l=xx*x_velocity(p_cell%l_state)+yy*y_velocity(p_cell%l_state)
    u_radial_r=xx*x_velocity(p_cell%r_state)+yy*y_velocity(p_cell%r_state)
    u_radial=xx*x_velocity(p_cell%or_state)+yy*y_velocity(p_cell%or_state)
    u_tangen_l=-yy*x_velocity(p_cell%l_state)+xx*y_velocity(p_cell%l_state)
    u_tangen_r=-yy*x_velocity(p_cell%r_state)+xx*y_velocity(p_cell%r_state)
    u_tangen=-yy*x_velocity(p_cell%or_state)+xx*y_velocity(p_cell%or_state)
    
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

! call scanning(ggd_cell)

do i=1,nodes_number
 if(ndd(i)%status.eq.'awake') then
  n_cell=ndd(i)%n_cell
  uu_show(n_cell%x_idx,n_cell%y_idx)=n_cell%or_state
 end if
end do

! sum=0.d0
! do i=nxll, nxll+nx-1
!  do j=nyll, nyll+ny-1
!   sum=sum+h*h*u_show(i,j)
!  end do
! end do

do i=nxll, nxll+nx-1
 do j=nyll,nyll+ny-1
  rho_show(i,j)=density(uu_show(i,j))
  u_show(i,j)=x_velocity(uu_show(i,j))
  v_show(i,j)=y_velocity(uu_show(i,j))
  if(ggd_cell(i,j)%region.eq.'smth') then
   p_show(i,j)=pressure(uu_show(i,j))
  end if
 end do
end do

do i=1,curves_number
 if(cvv(i)%status.eq.'awake') then
  nullify(temp)
   
  if(associated(cvv(i)%begin%next)) then
   temp=>cvv(i)%begin%next
   head_mark=temp%address%idx; head_switch=.true.
  end if
  do while(associated(temp).and.head_switch) 
   g_cell=temp%g_cell; p_cell=temp%p_cell
   ix=g_cell%x_idx; jy=g_cell%y_idx
   if(p_cell%wv_nb.eq.2) then
    p_show(ix,jy)=pressure(p_cell%r_state)
   end if
        
   temp=>temp%next
   if(associated(temp)) then
    head_switch=(temp%address%idx.ne.head_mark)
   else
    exit
   end if
  end do
 end if
end do

! do i=nxll, nxll+nx-1
!  do j=nyll,nyll+ny-1
!   xx=dfloat(i); yy=dfloat(j)
!   rr=dsqrt(xx*xx+yy*yy)
!   if(rr.gt.0.1d-6) then
!    xx=xx/rr; yy=yy/rr
!   else
!    xx=0.0d0; yy=0.0d0
!   end if
!   velo_radial(i,j)=u_show(i,j)*xx+v_show(i,j)*yy
!   velo_tangential(i,j)=-yy*u_show(i,j)+xx*v_show(i,j)     
!  end do
! end do

! write(*,*) ' Please input I0 and J0 '
! read(*,'(2i5)') i0, j0
! call scan_shw(1,i0,j0)
! call scan_shw(2,i0,j0)
! call scan_shw(3,i0,j0) 
! call scan_shw(4,i0,j0)

! call scan_shw(5,i0,j0) 

! call scan_shw(2,-6,7) 
! call scan_shw(3,-6,7)

! Compute vorticity.

! xxxx=0.0d0; yyyy=0.0d0
do i=nxll+1,nxll+nx-2
 do j=nyll+1,nyll+ny-2
  vorticity_show(i,j)=(v_show(i+1,j)-v_show(i-1,j))/h/2.0d0
  vorticity_show(i,j)=vorticity_show(i,j)-(u_show(i,j+1)-u_show(i,j-1))/h/2.0d0

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! Compute vorticities for Jacobs-like.
!   if(dfloat(i)*h.ge.1.0d0) then        
!    if(vorticity_show(i,j).gt.0.0d0) then
!     xxxx=xxxx+vorticity_show(i,j)*h*h
!    else
!     yyyy=yyyy+vorticity_show(i,j)*h*h
!    end if
!   end if
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! Compute vorticities for Air-He and Air-R22.
!   if(dfloat(j)*h.ge.0.0d0) then        
!    if(vorticity_show(i,j).gt.0.0d0) then
!     xxxx=xxxx+vorticity_show(i,j)*h*h
!    else
!     yyyy=yyyy+vorticity_show(i,j)*h*h
!    end if
!   end if
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!   vorticity_show(i,j)=dabs(vorticity_show(i,j))
 end do
end do
 
do i=nxll+1,nxll+nx-2
 vorticity_show(i,nyll)=vorticity_show(i,nyll+1)
 vorticity_show(i,nyll+ny-1)=vorticity_show(i,nyll+ny-2)
end do

do j=nyll+1,nyll+ny-2
 vorticity_show(nxll,j)=vorticity_show(nxll+1,j)
 vorticity_show(nxll+nx-1,j)=vorticity_show(nxll+nx-2,j)
end do

vorticity_show(nxll,nyll)=vorticity_show(nxll+1,nyll+1)
vorticity_show(nxll+nx-1,nyll)=vorticity_show(nxll+nx-2,nyll+1)
vorticity_show(nxll,nyll+ny-1)=vorticity_show(nxll+1,nyll+ny-2)
vorticity_show(nxll+nx-1,nyll+ny-1)=vorticity_show(nxll+nx-2,nyll+ny-2)

! aaa=0.0d0; bbb=0.0d0
! do i=nxll+1, nxll+nx-2
!  do j=nyll+1, nyll+ny-2
!   if(vorticity_show(i,j).lt.aaa) aaa=vorticity_show(i,j)
!   if(vorticity_show(i,j).gt.bbb) bbb=vorticity_show(i,j)
!  end do
! end do

!! Compute Mach number.
!do i=nxll,nxll+nx-1
! do j=nyll,nyll+ny-1
!  call compute_sound_speed(i,j,bbb)
!  aaa=u_show(i,j)-0.4d0 ! Tentatively for airHe shock-bubble.
!  mach_number(i,j)=aaa/bbb
! end do
!end do
 
!! Compute acoustic_impedance.
!do i=nxll,nxll+nx-1
! do j=nyll,nyll+ny-1
!  call compute_sound_speed(i,j,bbb)
!  acoustic_impedance(i,j)=bbb*rho_show(i,j)
! end do
!end do

! call scan_shw(5,i0,j0)

open(1, file='d:\workf90_1\show\shown.dat')
write(1,'(f10.4)') x_width, y_width
write(1,'(i5)') nxll, nyll, nx, ny
close(1)

open(2, file='d:\workf90_1\show\showrho.dat')
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  write(2,'(f12.6)') rho_show(i,j)
 end do
end do  
close(2)
! Ouput for show of solution.

open(2, file='d:\workf90_1\show\showu.dat')
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
!  write(2,'(f12.6)') uu_show(i,j)%value(2)
  write(2,'(f12.6)') u_show(i,j)
 end do
end do  
close(2)

open(2, file='d:\workf90_1\show\showv.dat')
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
!  write(2,'(f12.6)') uu_show(i,j)%value(3)
  write(2,'(f12.6)') v_show(i,j)
 end do
end do  
close(2)

open(2, file='d:\workf90_1\show\showp.dat')
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
!  write(2,'(f12.6)') uu_show(i,j)%value(4)
  write(2,'(f12.6)') p_show(i,j)
 end do
end do  
close(2)

open(6, file='d:\workf90_1\show\showvort.dat')
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
!  write(2,'(f12.6)') uu_show(i,j)%value(4)
  write(6,'(f12.6)') vorticity_show(i,j)
 end do
end do  
close(6)

open(6, file='d:\workf90_1\show\showmach.dat')
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
!  write(2,'(f12.6)') uu_show(i,j)%value(4)
  write(6,'(f12.6)') mach_number(i,j)
 end do
end do  
close(6)

open(6, file='d:\workf90_1\show\showacoustic.dat')
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
!  write(2,'(f12.6)') uu_show(i,j)%value(4)
  write(6,'(f12.6)') acoustic_impedance(i,j)
 end do
end do  
close(6)

! output the show of solution.

! open(2, file='d:\workf90_1\show\velo_r.dat')
! do i=nxll,nxll+nx-1
!  do j=nyll,nyll+ny-1
!   write(2,'(f12.6)') velo_radial(i,j)
!  end do
! end do  
! close(2)

! open(2, file='d:\workf90_1\show\velo_t.dat')
! do i=nxll,nxll+nx-1
!  do j=nyll,nyll+ny-1
!   write(2,'(f12.6)') velo_tangential(i,j)
!  end do
! end do  
! close(2)

!call check_list(cvv(1),'down  ')

call curve_files(file)
m=0
do lc=1,curves_number
 if(cvv(lc)%status.eq.'awake') then
  m=m+1; ll=0 
  temp=>cvv(lc)%begin%next
  head_mark=temp%address%idx; head_switch=.true.
  if(associated(temp)) then
   g_cell=temp%g_cell

   select case(g_cell%edge(1))
    case(1)
     x(0)=(g_cell%x_idx+g_cell%dis_posi(1))*h
     y(0)=(g_cell%y_idx-0.5d0)*h
   	case(2)
     x(0)=(g_cell%x_idx+0.5d0)*h 
     y(0)=(g_cell%y_idx+g_cell%dis_posi(1))*h
    case(3)
     x(0)=(g_cell%x_idx+g_cell%dis_posi(1))*h
     y(0)=(g_cell%y_idx+0.5d0)*h
   	case(4)
     x(0)=(g_cell%x_idx-0.5d0)*h
     y(0)=(g_cell%y_idx+g_cell%dis_posi(1))*h
   end select

   do while(associated(temp).and.head_switch)
    g_cell=temp%g_cell; ll=ll+1
    select case(g_cell%edge(2))
     case(1)
      x(ll)=(g_cell%x_idx+g_cell%dis_posi(2))*h
      y(ll)=(g_cell%y_idx-0.5d0)*h
     case(2)
      x(ll)=(g_cell%x_idx+0.5d0)*h
      y(ll)=(g_cell%y_idx+g_cell%dis_posi(2))*h
     case(3)
      x(ll)=(g_cell%x_idx+g_cell%dis_posi(2))*h
      y(ll)=(g_cell%y_idx+0.5d0)*h
     case(4)
      x(ll)=(g_cell%x_idx-0.5d0)*h 
      y(ll)=(g_cell%y_idx+g_cell%dis_posi(2))*h
    end select
    temp=>temp%next
	if(associated(temp))head_switch=(temp%address%idx.ne.head_mark)
   end do

  end if
 
  open(4, file=file(m))
  do l=0,ll
   write(4,'(2f12.6)') x(l), y(l)
  end do
  close(4)
 end if
end do

open(5, file='d:\workf90_1\show\showm.dat')
write(5,'(i5)') m
close(5)
! Output for the show of discontinuity curves.

deallocate(uu); deallocate(cvv)
deallocate(ggd_cell); deallocate(ndd)
deallocate(uu_show)
deallocate(rho_show); deallocate(u_show); deallocate(v_show)
deallocate(p_show)
deallocate(velo_radial); deallocate(velo_tangential)
deallocate(vorticity_show)
deallocate(mach_number)
deallocate(acoustic_impedance)

end subroutine show_out


subroutine curve_files(file)
implicit none
character*38 :: file(10)

file(1)='d:\workf90_1\show\showdc1.dat'
file(2)='d:\workf90_1\show\showdc2.dat'
file(3)='d:\workf90_1\show\showdc3.dat'
file(4)='d:\workf90_1\show\showdc4.dat'
file(5)='d:\workf90_1\show\showdc5.dat'
file(6)='d:\workf90_1\show\showdc6.dat'
file(7)='d:\workf90_1\show\showdc7.dat'
file(8)='d:\workf90_1\show\showdc8.dat'
file(9)='d:\workf90_1\show\showdc9.dat'
file(10)='d:\workf90_1\show\showdc10.dat'
end subroutine curve_files


subroutine scan_shw(quantity,i0,j0)

implicit none
integer, intent(in) :: quantity,i0,j0

real(8), dimension(-3:3,-3:3) :: uuu
integer :: ii, jj

! i0=50; j0=72

select case(quantity)
 case(1)
  print*, 'Density ='
  uuu=rho_show(i0-3:i0+3,j0-3:j0+3)
  print*, ''
 case(2)
  print*, 'X-velocity ='
  uuu=u_show(i0-3:i0+3,j0-3:j0+3)
  print*, ''
 case(3)
  print*, 'Y-velocity ='
  uuu=v_show(i0-3:i0+3,j0-3:j0+3)
  print*, ''
 case(4)
  print*, 'Pressure ='
  uuu=p_show(i0-3:i0+3,j0-3:j0+3)
  print*, ''
 case(5)
  print*, 'Vorticity ='
  uuu=vorticity_show(i0-3:i0+3,j0-3:j0+3)
  print*, ''
 case default; call error_message
end select
       
do jj=3,-3,-1
 write(*,'(7f10.6)') (uuu(ii,jj), ii=-3,3)
end do

end subroutine scan_shw


subroutine compute_sound_speed(i,j,sd)

implicit none

integer, intent(in) :: i, j
real(8), intent(out) :: sd

real(8) :: sdl, sdr, areal, arear
type(critical_cell), pointer :: temp_g


select case(ggd_cell(i,j)%region)

 case('smth  '); sd=sound_speed(uu(i,j))

 case('crit  ')
  call visit(ggd_cell(i,j)%ccpt(1)%address,temp_g)
  sdl=sound_speed(temp_g%p_cell%l_state)
  sdr=sound_speed(temp_g%p_cell%r_state)
  areal=side_area(temp_g%g_cell,'left  ')
  arear=side_area(temp_g%g_cell,'right ')
  sd=areal*sdl+arear*sdr

 case default; call error_message

end select

end subroutine compute_sound_speed

end module output_show