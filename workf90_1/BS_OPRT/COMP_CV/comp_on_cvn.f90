module computation_on_curves
! This module carries out the computation on auxiliary discontinuity 
! curves.

!use wind_directions
! 'windn.f90'

!use corner_check
! 'ch_crnsn.f90'

use compute_cell_average
! 'comp_can.f90'

use recover_mid_positions
! 'rcv_mpn.f90'

use compute_discontinuity_positions
! 'comp_dpn.f90'
use compute_boundary_positions
! 'bd_posin.f90'

use restriction_of_positions
! 'restr_dp.f90'

use half_step_computation
! 'half_stp.f90'

use compute_node_ca
! 'comp_ndn.f90'

use output_middle_positions
! 'out_mdl.f90'

use boundary_conditions
! 'boundary_conditions.f90'

use numerical_regularization
! 'num_regular.f90'

use update_left_n_right_states
! 'updat_lr.f90'

!use manipulation_1
!! 'manipu_1.f90'

! use manipulation_3
!! 'manipu_3.f90'
 
implicit none

public  compute_cell_averages_on_curves
private recovery, compute_positions, restrict_positions, half_step, &
        compute_cvcas, boundary_positions, comput_ndca !, manipulate_1 


contains


subroutine compute_cell_averages_on_curves

implicit none

! type(auxi_crit_cell), pointer :: tempttt, tempsss, temp000
! integer :: llllh

allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(acvv(cvn))
allocate(ndd(ndn))
call give_acv; call give_nd; call give_gcl

call give_sth(0)

!call check_corners
!call compute_winds

! call check_list_ac(1,'up    ')

call data_shift

! call check_list_ac(1,'down  ')

! call check_partner_conserv_ac(1)

!call manipulate_1

! call check_list_ac(1,'down  ') 

call recovery

! call out_middle_posi(1)
! call check_list_ac(1,'down  ')
 
! pause

call compute_positions

! call check_list_ac(1,'down  ')

call boundary_condition3
call boundary_condition2

! call check_list_ac(1,'down  ')

call restrict_positions(0)

!call numerical_regularize

! call check_list_ac(1,'down  ')

call half_step

! call check_list_ac(1,'down  ')

call boundary_condition3
call boundary_condition2

! call check_list_ac(1,'down  ')

call restrict_positions(1)

! call check_positions_1(1)

! call check_partner_conserv_ac(1)

! call manipulate_2

call compute_cvcas

! call check_list_ac(1,'down  ')

call data_shift

! call check_partner_conserv_ac(1)

call give_sth(1)

! call check_partner_conserv_ac(1)

call boundary_condition2

! call manipulate_1

call recovery

! call manipulate_3

! call out_middle_posi(1)

! call check_list_ac(1,'down  ')

! call check_partner_conserv_ac(1)

call compute_positions

call boundary_condition3
call boundary_condition2

! call check_list_ac(1,'down  ')

!call boundary_positions

call restrict_positions(1) 

!call boundary_condition3
!call boundary_condition2

! call check_list_ac(1,'down  ')

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!call update_left_n_right
! Update the left and right states in call critical cells to let
! information associated with other characteristic fields propagate
! across the tracked discontinuity.

!Temporary stop.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! call check_partner_conserv_ac(1)

! call check_list_ac(1,'down  ')

call boundary_condition2

! temp000=>acvv(1)%begin%next 
! if(temp000%g_cell%y_idx.ne.0) then
!  temp000=>temp000%previous
!  if(temp000%g_cell%y_idx.ne.0) then
!   temp000=>acvv(1)%eend%previous
!   if(temp000%g_cell%y_idx.ne.0) call error_message
!  end if
! end if
! temp000%g_cell%dis_posi(1)=temp000%g_cell%dis_posi(2)
! tempttt=>temp000%next; tempsss=>temp000%previous
! do llllh=1, 7
!  tempsss%g_cell%dis_posi(1)=tempttt%g_cell%dis_posi(2)
!  tempsss%g_cell%dis_posi(2)=tempttt%g_cell%dis_Posi(1)
!  tempttt=>tempttt%next; tempsss=>tempsss%previous
! end do

! call check_list_ac(1,'down  ')

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
call numerical_regularize

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! call check_list_ac(1,'down  ')

! call check_partner_conserv_ac(1)

! call manipulate_2

call comput_ndca


call get_acv; call get_nd
deallocate(ggd_cell); deallocate(acvv); deallocate(ndd); deallocate(uu)


contains


subroutine check_positions_1(cv_nmb)

implicit none
integer, intent(in) :: cv_nmb

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(comp_info) :: c_cell
integer :: head_mark, end_mark
logical :: head_switch

temp=>acvv(cv_nmb)%begin%next
head_mark=temp%address%idx; head_switch=.true.
head_mark=temp%address%idx; head_switch=.true.
boundary_end=>acvv(cv_nmb)%eend_boundary%previous
if(acvv(cv_nmb)%eend%previous%address.eq.boundary_end%address) then
 end_mark=error_index
else
 end_mark=acvv(cv_nmb)%eend%previous%next%address%idx
end if

do while(associated(temp).and.head_switch)
 g_cell=temp%g_cell; c_cell=temp%c_cell
 g_cell%dis_posi=0.5d0*(c_cell%tdp+g_cell%dis_posi)

 temp%g_cell=g_cell

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if

end do

end subroutine check_positions_1


end subroutine compute_cell_averages_on_curves


subroutine data_shift
! The left and right states in auxiliary critical cells are shifted between the initial
! and final time levels, the later ones are stored in 'c_cell%l_state(0)'s and 
! 'c_cell%r_state(0)'s.

implicit none

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell
type(phy_info) :: p_cell
type(comp_info) :: c_cell
integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake ') then
  call shift(acvv(i))
 end if
end do


contains


subroutine shift(acv)

implicit none

type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch
type(state) :: swap

nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
 boundary_end=>acv%eend_boundary%previous
 if(acv%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acv%eend%previous%next%address%idx
 end if
end if

do while(associated(temp).and.head_switch)

 address=temp%address
 g_cell=temp%g_cell; p_cell=temp%p_cell; c_cell=temp%c_cell
 
 swap=p_cell%l_state
 p_cell%l_state=c_cell%left_room_for_swap
 c_cell%left_room_for_swap=swap
 swap=p_cell%r_state
 p_cell%r_state=c_cell%right_room_for_swap
 c_cell%right_room_for_swap=swap

 temp%p_cell=p_cell; temp%c_cell=c_cell

 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
end do

end subroutine shift


end subroutine  data_shift


end module computation_on_curves