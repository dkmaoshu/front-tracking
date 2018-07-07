program main

!use compute_growth_rate
!new

!use save_growth_rate
!new

use input_setup
! 'input?.f90'

use out_intermediate
! 'output?.f90'

use output_show
! 'output_for_show.f90

use auxiliary_curves_production
! 'auxi_cvp.f90'

use computation_in_smooth
! 'comp_in_smooth.f90'

use computation_on_curves
! 'comp_on_ca.f90'

use node_move_n_split
! 'nd_mvsp.f90'

use reset_discvs
! 'reset.f90'

use move_all
! 'mv_all.f90'

use mesh_ratio
! 'mesh_rt.f90'	

use compress_linear_dis
! 'comlids.f90.f90'

use output_show_HR
! 'sg_sts_sc_HR.f90'

use revers_discontinuity_curves
! 'rev_cvs.f90'

use compute_circulation_HS
! 'comp_circ_HS.f90'
 
!  real(8) circulation


!real(8) :: final_time
integer :: time_step=1, maximum_step

call input_and_setup

print*, 'Input the maximum time step'
read(*,'(i5)') maximum_step
print*, 'Input the final time'
read(*,'(f18.12)') final_time

! call revers_curves

! open(16,file='d:\workf90_1\show\circulation_HS.dat')

do while(current_time.lt.final_time.and.time_step.le.maximum_step) 
  
 total_step=total_step+1
  
 print*, ' time_step =', time_step, '   current_time =', current_time
 
 call mesh_ratio_compute

!  r=0.0d0

! call get_growth_rate(hal_amplitude,growth_rate)
!new
 
! call save_gr(time_step,amp,grt)
!new
 
 call production(auxiliary_curve_type)
 call compute_smooth
 call compute_cell_averages_on_curves
 call reset_discontinuity_curves  
 call move_n_splits
 call compress_linear_discontinuities

 call change_auxi_type
 
 current_time=current_time+r*h
 time_step=time_step+1

!  if(imod(time_step,20).eq.0) then
!   call compute_circulation(circulation)  
!   write(16,'(2f18.12)') current_time*dsqrt(10.0d0)*100.0d0, circulation
!  end if
  
 if(boundary_type_1.eq.'periodic    ') call move_region(.false.)

 if(mod(total_step,2).eq.0) call revers_curves

!  call show_out
!  pause
   
!  if(mod(time_step,10).eq.0) then
!   call output_m
!   call show_out
! end if
   
end do

! call compute_circulation(circulation)  
! write(16,'(2f18.12)') current_time*dsqrt(10.0d0)*100.0d0, circulation
 
! close(16)

if(boundary_type_1.eq.'periodic    ') call move_region(.true.)

call output_m
call show_out
!call show_out_hr


!open(4,file='d:\workf90_1\output\gr.dat')
!write(4,'(f18.12)') amp
!write(4,'(f18.12)') grt
!close(4)

!print*, amp
!print*, grt

end program main