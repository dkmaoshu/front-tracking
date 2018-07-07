module reset_discvs
! This module is for reseting the auxiliary discontinuity curves
! and transform them back to ordinary discontinuity curves.

!use update_left_n_right_states
!! 'updat_lr.f90'

use produce_local_smooth_4_nodes
! 'pro_smnd.f90'

use local_smooth_fix
! 'lcl_fix.f90'

use data_preparation
! 'data_pre.f90'

use update_curves
! 'upd_cvs.f90'

use fix_reversed_curves
! 'fix_recv.f90'

use form_discurves
! 'form_dsv.f90'

use clean_up
! 'clean_up.f90'

use disjoints_fix
! 'dis_jnt.f90'

use smoothen_curves
! 'smooth_curves.f90' 

use packed_separation
! 'packed_separate.f90'

!use manipulation_2
!! 'manipu_2n.f90'

!use manipulation_3
!! 'manipu_3.f90'

use manipulation_4
! 'manipu_3.f90'

use boundary_conditions
! 'boundary_conditions.f90'

use numerical_regular_partners
! 'num_regular.stack1.f90'

use numerical_NS_diffusion
! 'NS_diffusion.f90'

implicit none

!type(crit_info) :: crit

public  reset_discontinuity_curves
private uu, acvv, ggd_cell, cvv, update_acvs, ndd, form_new, &
        fix_sides, smnd_fix, change_order_of_curve_reset, &
		clean_up_discontinuities
        

contains


subroutine reset_discontinuity_curves

implicit none
integer :: i
 
! type(grid_cell) :: ggd_smpl
!type(crit_info) :: ccc
		
allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(acvv(curves_number)) 
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(curves_number))
allocate(ndd(ndn))

call give_sth(1)
call give_acv
call give_gcl
call give_nd

!ccc=crit

! call scanning(uu)
! call check_ring(1,'counter ')
! call check_list_ac(1,'down  ')

! call check_partner_conserv_ac(1)

! ggd_smpl=ggd_cell(-1,3)

!call update_left_n_right
!! Update the left and right states in call critical cells to let
!! information associated with other characteristic fields propagate
!! across the tracked discontinuity.

! call check_partner_conserv_ac(1)

call boundary_condition2

do i=1, ndn
 if(ndd(i)%status.eq.'awake ') then
  call prlcsmnd(i)
 end if
end do

! call check_ring('counter ')

! call check_list_ac(1,'down  ')

do i=1, ndn
 if(ndd(i)%status.eq.'awake ') then
  call smnd_fix(i)
 end if
end do

! call check_ring('counter ')
! call check_list_ac(1,'down  ')

call data_preparing

call boundary_condition2

call make_partner_memory

! call check_list_ac(1,'down  ')

call update_acvs
! Update the auxiliary discontinuity curves by moveing and 
! splitting auxiliary critical cells if necessary.

! call check_list_ac(1,'up    ')

call numerical_regularize_p
! Or_states of splitted discontinuity cells are recomputed here. The computation of or_states in "update_acvs" is not
! quite stable.

! call check_list_ac(1,'down  ')

call manipulate_4

! call check_list_ac(1,'down  ')

call fix_reversed_cv

! call scanning(ggd_cell)
! call check_list_ac(1,'down  ')
! pause

call form_new
! Transform the auxiliary discontinuity curves back to the 
! ordinary discontinuity curves.

! call check_list_c(1,'down  ')

!call smoothen_cv

call clean_up_discontinuities

! call check_list_c(1,'down  ')

call num_NS_diffusion
! Numerical diffusion in Navier-Stokes form.

! call check_list_c(1,'down  ')

! Clean again after the numerical Navier-Stokes diffusion.
call clean_up_discontinuities

! call check_list_c(1,'up    ')

call fix_all_disjoints

call fix_sides

! call check_list_c(1,'down  ')

! call manipulate_2
! call check_list_c(1,'down  ')

call packed_separate

!call manipulate_3

! call check_list_c(1,'down  ')
! pause

call curves_type_change('regula')

call change_order_of_curve_reset

! ggd_smpl=ggd_cell(3,3)

! call scanning(uu)

call get_sth(0)
call get_cv
call get_gcl
call get_nd

deallocate(uu); deallocate(acvv); deallocate(ggd_cell)
deallocate(cvv); deallocate(ndd)

end subroutine reset_discontinuity_curves


subroutine fix_sides
! This subroutine is written for fixing side-problems in the grid-map at the end of 
! reset step.

implicit none

type(critical_cell), pointer :: temp
type(geo_info) :: g_cell
type(adss_info) :: address
integer :: head_mark, i0, j0, l, i
logical :: head_switch

do l=1, cvn
 if(cvv(l)%status.eq.'awake ') then
  if(associated(cvv(l)%begin%next)) then
   temp=>cvv(l)%begin%next
   head_mark=temp%address%idx; head_switch=.true.
   do while(associated(temp).and.head_switch)
    g_cell=temp%g_cell; address=temp%address
    i0=g_cell%x_idx; j0=g_cell%y_idx
    do i=1,4
     if(ggd_cell(i0,j0)%ccpt(i)%address.eq.address) then
      if(ggd_cell(i0,j0)%ccpt(i)%side.ne.g_cell%point(i)) then
       print*, 'A side-fix is being done', l, i0, j0, i, ggd_cell(i0,j0)%ccpt(i)%side
       ggd_cell(i0,j0)%ccpt(i)%side=g_cell%point(i)
      end if
     end if
    end do

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

end subroutine fix_sides
 
end module reset_discvs