program Peng_Zabusky_Zhang_problem

! Gaozhu Peng, Norman J. Zabusky and Shuang Zhangs' problem, 
! Physics of Fluids, Vol.15 No.12, Dec. (2003).

use form_solution
! 'form_sl.f90'

integer :: curves, ix, jy

! type(grid_cell) :: ggd_c

print*, 'Please input the number of grid cells, NN ='
read(*,'(i5)') nn

print*, 'Please choose a case, CASES ='
read(*,'(i5)') cases

 !nxll=0; nyll=-nn; nx=nn+1; ny=2*nn+1
 !nxll=0; nyll=-8*nn; nx=nn+1; ny=16*nn+1  !(Air-SF6)
 !nxll=0; nyll=-0.5*(30/4)*nn; nx=nn+1; ny=(30/4)*nn+1   !(Air-He)
 
 
 !nxll=0; nyll=-9*nn; nx=nn+1; ny=27*nn/2+1  ! Modified Peng et al.(Air-SF6).

 nxll=0; nyll=-nn; nx=nn+1; ny=2*nn+1  ! Modified Peng et al.(Air-SF6) with nonreflecting boundaries at the two ends.

 !! SECOND PENG!!!
 !nxll=0; nyll=-7*nn/2; nx=nn+1; ny=7*nn+1

 ! x_width=1.0d0; y_width=2.0d0      ! Mushroom case
 ! x_width=3.75d0; y_width=22.5d0    !(Air-SF6, Holmes)
 ! x_width=3.75d0; y_width=60.0d0    !(Air-SF6, Ours)
 ! x_width=4.0d0; y_width=30.0d0     !(Air-He)
 x_width=2.0d0; y_width=4.0d0      ! Modified Peng et al (Air-SF6)


call set_grid
call set_sth
call set_gcl


allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(curves_number))
allocate(ndd(nodes_number))

current_time=0.d0
do ix=nxll-3,nxll+nx+2; do jy=nyll-3,nyll+ny+2 
 uu(ix,jy)=error_data
end do; end do

! First, set the solution in smooth region.
call set_in_smooth(cases)

! Then, set the contact discontinuity curve.
call set_discontinuity_curve

! call check_list_c(1,'down  ') 

do ii=2, cvn
 cvv(ii)%status='asleep'
end do

auxiliary_curve_type='xx '; order_of_curve_reset='inverse_way '  

! Initially, the auxiliary discontinuity curves will be built in
! 'xx'-type.

do ii=1, ndn
 ndd(ii)%status='asleep'
end do

! call scanning(ggd_cell)
! call check_list_c(1,'down  ')

boundary_type_1='Jacobs_PengZ'; boundary_type_2='nonreflect  '
total_step=0

call get_sth(0)
call get_cv
call get_gcl
call get_nd

deallocate(uu); deallocate(cvv)
deallocate(ggd_cell); deallocate(ndd)

call output_m
print*, ''

end program Peng_Zabusky_Zhang_problem

