program HS_He_air_1d

! Haas-Sturtevant shock-bubble interaction experiments, air-helium case, one-dimensional case. 

use form_solution_1d
! 'form_sl.f90'

integer :: curves, ix, jy

! type(grid_cell) :: ggd_c

print*, 'Please input the number of grid cells, NN ='
read(*,'(i5)') nn

print*, 'Please choose the radius, RADIUS ='
read(*,'(f18.12)') radius

print*, 'Please choose the case, CASES ='
read(*,'(i5)') cases

!nxll=-5*nn; nyll=-nn; nx=10*nn+1; ny=2*nn+1
!x_width=4.450d0; y_width=0.89d0

nxll=-3*nn; nyll=-5; nx=6nn+1; ny=11
x_width=2.67d0

!nxll=-10*nn; nyll=-nn; nx=20*nn+1; ny=2*nn+1
!x_width=8.90d0; y_width=0.89d0

call set_grid; y_width=h*10.0d0
call set_sth
call set_gcl

allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(curves_number))
allocate(ndd(nodes_number))

current_time=0.d0
do ix=nxll,nxll+nx-1; do jy=nyll,nyll+ny-1 
 uu(ix,jy)=error_data
end do; end do

! First, set the solution in smooth region.
call set_in_smooth

! Then, set the contact discontinuity curve.
!curves=1

call set_discontinuity_curve

! call check_list_c(2,'down  ')

!do curves=1,1
! call set_slops_on_curve(curves)
!end do

!do curves=1,1
! call compute_or_states(curves)
!end do

do ii=3, cvn
 cvv(ii)%status='asleep'
end do

auxiliary_curve_type='xx '; order_of_curve_reset='right_way   '

! Initially, the auxiliary discontinuity curves will be built in
! 'xx'-type.

do ii=1, ndn
 ndd(ii)%status='asleep'
end do

! call scanning(ggd_cell)
! call check_list_c(1,'down  ')

boundary_type='Haas_Stur_1d'
total_step=0

call get_sth(0)
call get_cv
call get_gcl
call get_nd

deallocate(uu); deallocate(cvv)
deallocate(ggd_cell); deallocate(ndd)

call output_m
print*, ''


end program HS_He_air_1d


