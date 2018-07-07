module mesh_ratio
! This module is for the computation of mesh ratio.

use solution
! 'solution.f90'

use numerical_fluxes, max_chsp=>maximum_characteristic_speeds
! 'flux_weno.f90'

real(8) :: global_courant_number
type(state) :: ch_sp

public  mesh_ration_compute
private global_courant_number


contains


subroutine mesh_ratio_compute
! Compute the mesh ratio.

implicit none

real(8) :: time_difference

allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(curves_number))
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
call give_sth(0)
call give_cv
call give_gcl

max_chsp=0.0d0
call courant_number_in_smooth
call courant_number_on_curves
if(global_courant_number.le.0.0d0) call error_message
r=rr/global_courant_number

time_difference=final_time-current_time
if(r*h.gt.time_difference) r=time_difference/h
max_chsp=1.1d0*max_chsp; max_chsp%value(3)=max_chsp%value(2)

deallocate(uu)
deallocate(cvv)
deallocate(ggd_cell)


contains 


subroutine courant_number_in_smooth

implicit none

integer :: i, j, k
real(8) :: crn, x_posi, y_posi

global_courant_number=0.0d0

do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  if(ggd_cell(i,j)%region.eq.'smth') then
   x_posi=dfloat(i)*h; y_posi=dfloat(j)*h
   call courant_number(uu(i,j),x_posi, y_posi,ch_sp,crn)
   do k=1, 4 
    max_chsp%value(k)=dmax1(max_chsp%value(k),ch_sp%value(k))
   end do  
   global_courant_number=dmax1(crn,global_courant_number)
  end if
 end do
end do

end subroutine courant_number_in_smooth


subroutine courant_number_on_curves

implicit none

type(critical_cell), pointer :: temp
integer :: head_mark, k, l
integer :: head_switch
type(geo_info) :: g_cell
type(phy_info) :: p_cell
real(8) :: crn, x_posi, y_posi

do k=1,curves_number
 if(cvv(k)%status.eq.'awake') then
  if(associated(cvv(k)%begin%next)) then
   temp=>cvv(k)%begin%next
   head_mark=temp%address%idx; head_switch=.true.

   do while(associated(temp).and.head_switch)
    
    g_cell=temp%g_cell; p_cell=temp%p_cell
	x_posi=dfloat(g_cell%x_idx)*h; y_posi=dfloat(g_cell%y_idx)*h
    call courant_number(p_cell%l_state,x_posi,y_posi,ch_sp,crn)
    do l=1, 4 
     max_chsp%value(l)=dmax1(max_chsp%value(l),ch_sp%value(l))
    end do  
    global_courant_number=dmax1(crn,global_courant_number)
	call courant_number(p_cell%r_state,x_posi,y_posi,ch_sp,crn)
    do l=1, 4 
     max_chsp%value(l)=dmax1(max_chsp%value(l),ch_sp%value(l))
    end do  
    global_courant_number=dmax1(crn,global_courant_number)

    temp=>temp%next
    if(associated(temp)) then
     head_switch=(temp%address%idx.ne.head_mark)
    else
     exit
    end if
   end do
  end if
 end if
end do

end subroutine courant_number_on_curves


end subroutine mesh_ratio_compute


end module mesh_ratio
