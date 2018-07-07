program check_accuracy

use solution

use input_setup_f

use input_setup

implicit none

real(8), dimension(:,:,:), allocatable :: u_exact, u_numer 
real(8), dimension(:,:), allocatable :: error
character*6, dimension(:,:), allocatable :: side
integer :: nxll_exact, nyll_exact, nx_exact, ny_exact
integer :: nxll_numer, nyll_numer, nx_numer, ny_numer
integer :: i, j, ic, jc, head_mark, mod, cv_number
logical :: head_switch
type(critical_cell), pointer :: temp
character*17 :: orientation_exact, orientation
real(8) :: sum, xx, yy

! Input the exact solution on a very fine grid.
call input_and_setup_f

! call check_list_c(1,'up    ')

orientation_exact='counter-clockwise'

nxll_exact=nxll; nyll_exact=nyll
nx_exact=nx; ny_exact=ny
allocate(u_exact(nxll_exact:nxll_exact+nx_exact-1,nyll_exact:nyll_exact+ny_exact-1,2))
u_exact=error_data

do i=nxll_exact, nxll_exact+nx_exact-1
 do j=nyll_exact, nyll_exact+ny_exact-1
  if(ggd_cell(i,j)%region.eq.'smth  ') u_exact(i,j,1)=uu(i,j)%value(1)
 end do
end do

do cv_number=1, 2

 if(associated(cvv(cv_number)%begin)) then

  nullify(temp)
  temp=>cvv(cv_number)%begin%next
  head_mark=temp%address%idx; head_switch=.true.

  do while(associated(temp).and.head_switch) 
   i=temp%g_cell%x_idx; j=temp%g_cell%y_idx
   u_exact(i,j,1)=temp%p_cell%l_state%value(1)
   u_exact(i,j,2)=temp%p_cell%r_state%value(1)

   temp=>temp%next
   if(associated(temp)) then
    head_switch=(temp%address%idx.ne.head_mark)
   else
    exit
   end if
  end do

 end if

end do

do i=nxll_exact,nxll_exact+nx_exact-1
 do j=nyll_exact,nyll_exact+ny_exact-1
  if(u_exact(i,j,1).lt.0.9d0*error_data) call error_message
 end do
end do
   
do i=1,cvn
 call deletee(cvv(i))
end do

deallocate(uu); deallocate(uu1); deallocate(ggd_cell)

! Input the numerical solution on a coarse grid.
call input_and_setup

! call check_list_c(1,'down  ')

orientation='clockwise        '

nxll_numer=nxll; nyll_numer=nyll
nx_numer=nx; ny_numer=ny
allocate(u_numer(nxll_numer:nxll_numer+nx_numer-1,nyll_numer:nyll_numer+ny_numer-1,2))
u_numer=error_data
allocate(side(nxll_numer:nxll_numer+nx_numer-1,nyll_numer:nyll_numer+ny_numer-1))
side='one   '

do i=nxll_numer, nxll_numer+nx_numer-1
 do j=nyll_numer, nyll_numer+ny_numer-1
  if(ggd_cell(i,j)%region.eq.'smth  ') u_numer(i,j,1)=uu(i,j)%value(1)
 end do
end do

do cv_number=1, 2

 if(associated(cvv(cv_number)%begin)) then

  nullify(temp)
  temp=>cvv(cv_number)%begin%next
  head_mark=temp%address%idx; head_switch=.true.

  do while(associated(temp).and.head_switch) 
   i=temp%g_cell%x_idx; j=temp%g_cell%y_idx
   side(i,j)=side_of_point(temp%g_cell,(/0.0d0,0.0d0/))
!   select case(side)
!    case('left  '); u_numer(i,j)=temp%p_cell%l_state%value(1)
!    case('right '); u_numer(i,j)=temp%p_cell%r_state%value(1)
!    case default; call error_message
!   end select
   u_numer(i,j,1)=temp%p_cell%l_state%value(1)
   u_numer(i,j,2)=temp%p_cell%r_state%value(1)

   temp=>temp%next
   if(associated(temp)) then
    head_switch=(temp%address%idx.ne.head_mark)
   else
    exit
   end if
  end do

 end if

end do

do i=nxll_numer,nxll_numer+nx_numer-1
 do j=nyll_numer,nyll_numer+ny_numer-1
  if(u_numer(i,j,1).lt.0.9d0*error_data) call error_message
 end do
end do
   
do i=1,cvn
 call deletee(cvv(i))
end do

deallocate(uu); deallocate(ggd_cell)

mod=nxll_exact/nxll_numer
h=1.0d0/dfloat(nxll)

allocate(error(nxll_numer:nxll_numer+nx_numer-1,nyll_numer:nyll_numer+ny_numer-1))
error=error_data

sum=0.0d0
do i=nxll_numer,nxll_numer+nx_numer-1
 do j=nyll_numer, nyll_numer+ny_numer-1
  ic=i*mod; jc=j*mod
  select case(side(i,j))
   case('one   ') 
    if(u_exact(ic,jc,2).gt.0.9d0*error_data) call error_message
    error(i,j)=dabs(u_numer(i,j,1)-u_exact(ic,jc,1))  
    sum=sum+error(i,j)
   case('left  ')
    if(u_exact(ic,jc,2).lt.0.9d0*error_data) then
	 xx=dabs(u_numer(i,j,1)-u_exact(ic,jc,1))
	 yy=dabs(u_numer(i,j,2)-u_exact(ic,jc,1))
     error(i,j)=dmin1(xx,yy)
     sum=sum+error(i,j)
    else
     if(orientation.eq.orientation_exact) then
      error(i,j)=dabs(u_numer(i,j,1)-u_exact(ic,jc,1))
      sum=sum+error(i,j)
     else
      error(i,j)=dabs(u_numer(i,j,1)-u_exact(ic,jc,2))
      sum=sum+error(i,j)
     end if
    end if
   case('right')
    if(u_exact(ic,jc,2).lt.0.9d0*error_data) then
	 xx=dabs(u_numer(i,j,1)-u_exact(ic,jc,1))
	 yy=dabs(u_numer(i,j,2)-u_exact(ic,jc,1))
     error(i,j)=dmin1(xx,yy)
     sum=sum+error(i,j)
    else
     if(orientation.eq.orientation_exact) then
      error(i,j)=dabs(u_numer(i,j,2)-u_exact(ic,jc,2))
      sum=sum+error(i,j)
     else
      error(i,j)=dabs(u_numer(i,j,2)-u_exact(ic,jc,1))
      sum=sum+error(i,j)
     end if
    end if
  end select	   
 end do
end do

sum=sum*h*h

continue

 xx=3.117295672188016E-003
 yy=1.189684868696015E-003

 xx=(dlog10(xx)-dlog10(yy))/dlog10(2.0d0)


end program check_accuracy