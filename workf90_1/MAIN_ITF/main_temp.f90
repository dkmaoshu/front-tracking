program main

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
!! 'nd_mvsp.f90'

use reset_discvs
! 'reset.f90'

use move_all
!! 'mv_all.f90'

use mesh_ratio
! 'mesh_rt.f90'

use output_show_HR
! 'sh_sts_el_HR.f90'

real(8) :: final_time
integer :: time_step=1, maximum_step

call input_and_setup
print*, 'Input the maximum time step'
read(*,'(i5)') maximum_step
print*, 'Input the final time'
read(*,'(f18.12)') final_time

do while(current_time.lt.final_time.and.time_step.le.maximum_step) 

 print*, ' time_step =', time_step
 call mesh_ratio_compute
 call production(auxiliary_curve_type)
 call compute_smooth
 call compute_cell_averages_on_curves
 call reset_discontinuity_curves
 call move_n_splits

 call change_auxi_type
 current_time=current_time+r*h
 time_step=time_step+1

 if(boundary_type.eq.'periodic    ') call move_region

!  call show_out
!  pause

end do

call output_m
call show_out
!call show_out_hr

end program main