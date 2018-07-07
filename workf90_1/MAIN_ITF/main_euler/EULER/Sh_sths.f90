module show_smth

use physical_state
use euler_functions
implicit none


contains


subroutine show(u_sh)

implicit none

type(state), intent(in), dimension(-3:3,-3:3) :: u_sh
real(8), dimension(-3:3,-3:3) :: u_components
real(8), dimension(4) :: value
integer :: choice, i, j

print*, 'What do you want to see?'
read(*,'(i5)') choice

select case(choice)

 case(1,2,3,4)
  print*, 'So, you want to see the', choice, 'th components. Ok, here they are.'
  do i=-3,3
   do j=-3,3
    u_components(i,j)=u_sh(i,j)%value(choice)
   end do
  end do
 case(5)
  print*, 'So, you want to see the x-velocity. Ok, here they are.'
  do i=-3,3
   do j=-3,3
    value=u_sh(i,j)%value
    if(value(1).gt.0.9d0*error_data) then
     u_components(i,j)=uf(value)
    else
     u_components(i,j)=error_data
    end if	 
   end do
  end do
 case(6)
  print*, 'So, you want to see the y-velocity. Ok, here they are.'
  do i=-3,3
   do j=-3,3
    value=u_sh(i,j)%value
    if(value(1).gt.0.9d0*error_data) then
     u_components(i,j)=vf(value)
    else
     u_components(i,j)=error_data
    end if
   end do
  end do
 case(7)
  print*, 'So, you want to see the pressure. Ok, here they are.'
  do i=-3,3 
   do j=-3,3
    value=u_sh(i,j)%value
    if(value(1).gt.0.9d0*error_data) then
     u_components(i,j)=pf(value,u_sh(i,j)%gamma)
    else
     u_components(i,j)=error_data
    end if
   end do
  end do
end select

call project(u_components)

end subroutine show


subroutine project(u_components)

implicit none
real(8), intent(in), dimension(-3:3,-3:3) :: u_components
integer :: i, j

do j=3,-3,-1
 write(*,'(7f10.4)') (u_components(i,j),i=-3,3)
end do

end subroutine project


end module show_smth