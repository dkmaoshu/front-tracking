module form_discurves
! This module is for the formation of new discontinuity curves
! using information in auxiliary discontinuity curves and the 
! delete of the auxiliary ones after the formation.

use solu_comput
! 'solu_comp.f90'

use discontinuity_curves
! 'discontinuity_curves.f90'

implicit none

public form_new
private form_dc, init


contains


subroutine form_new
! Form all discontinuity curves and delete all auxiliary 
! discontinuity curves.
implicit none

integer :: i
 
! type(auxi_crit_cell), pointer :: temppx

do i=1,cvn
 if(acvv(i)%status.eq.'awake'.or.acvv(i)%status.eq.'yawn  ') then
  call init(acvv(i))
 end if
end do

! call check_list_ac(1,'down  ')

! Form discontinuity curves.
do i=1,cvn
 cvv(i)%status=acvv(i)%status
 if(acvv(i)%status.eq.'awake '.or.acvv(i)%status.eq.'yawn  ') then
  cvv(i)%cv_type=acvv(i)%cv_type
  cvv(i)%begin_end=acvv(i)%begin_end
  cvv(i)%end_end=acvv(i)%end_end
  cvv(i)%total=0; cvv(i)%wave=acvv(i)%wave
  call form_dc(i)
 end if
end do

! temppx=>acvv(2)%eend_boundary%previous
! call check_list_ac(2,'up    ')

! Delete auxiliary discontinuity curves.
do i=1,cvn
 call deletee(acvv(i))
end do

end subroutine form_new


subroutine init(acvv)

type(auxi_discv), intent(in) :: acvv

type(auxi_crit_cell), pointer :: tempa
integer :: head_mark
logical :: head_switch

nullify(tempa)
if(associated(acvv%begin_boundary%next)) then
 tempa=>acvv%begin%next
 head_mark=tempa%address%idx; head_switch=.true.
end if

do while(associated(tempa).and.head_switch) 

 call clean_up_address(tempa%transfered)

 tempa=>tempa%next
 if(associated(tempa)) then
  head_switch=(tempa%address%idx.ne.head_mark)
 else
  exit
 end if
end do

end subroutine init


subroutine form_dc(cv_nmb)
! Form a discontinuity curve 'cvv' from an auxiliary discontinuity
! curve 'acvv'.

integer, intent(in) :: cv_nmb

type(critical_cell), pointer :: temp, current, temp_g
type(auxi_crit_cell), pointer :: tempa, tempa_g, boundary_end
type(adss_info) :: adss_o, adss_n, adss_oo, adss_nn
integer :: head_mark, end_mark, i0, j0, i
logical :: head_switch
type(geo_info) :: g_cell, gn_cell
character*3, dimension(:,:,:), allocatable :: updated

allocate(updated(nxll:nxll+nx-1,nyll:nyll+ny-1,4))

updated='no '
call creat_cv(cv_nmb); if(acvv(cv_nmb)%status.eq.'yawn  ') return

current=>cvv(cv_nmb)%begin
nullify(tempa)
if(associated(acvv(cv_nmb)%begin%next)) then
 tempa=>acvv(cv_nmb)%begin%next
 head_mark=tempa%address%idx; head_switch=.true.
 boundary_end=>acvv(cv_nmb )%eend_boundary%previous
 if(acvv(cv_nmb)%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acvv(cv_nmb)%eend%previous%next%address%idx
 end if
end if

do while(associated(tempa).and.head_switch) 

 allocate(temp)
 temp%address%cv_nmb=tempa%address%cv_nmb
 g_cell=tempa%g_cell
 temp%g_cell=tempa%g_cell; temp%p_cell=tempa%p_cell
 temp%l_stk=tempa%l_stk; temp%r_stk=tempa%r_stk
 temp%tangled=tempa%tangled

 temp%g_cell%mdis_posi=error_data

!  call check_list(cvv,'down  ')

 call insert(current,temp,'after ')

!  call visit(adss_info(1,34),tt); gg=tt%next%g_cell

!  call check_list(cvv,'down  ')

 tempa%transfered=temp%address

! The following operation updates the spatial neighboring
! relations.
 adss_o=tempa%address; adss_n=temp%address

 adss_oo=tempa%l_stk
 if(adss_oo%cv_nmb.gt.0) then
  call visit(adss_oo,tempa_g)
  gn_cell=tempa_g%g_cell
  adss_nn=tempa_g%transfered
  if(adss_nn%cv_nmb.gt.0) then
   temp%l_stk=adss_nn
   call visit(adss_nn,temp_g)
   if(tempa_g%l_stk.eq.adss_o) temp_g%l_stk=adss_n
   if(tempa_g%r_stk.eq.adss_o) temp_g%r_stk=adss_n
  end if
 end if

 adss_oo=tempa%r_stk
 if(adss_oo%cv_nmb.gt.0) then
  call visit(adss_oo,tempa_g)
  gn_cell=tempa_g%g_cell
  adss_nn=tempa_g%transfered
  if(adss_nn%cv_nmb.gt.0) then
   temp%r_stk=adss_nn
   call visit(adss_nn,temp_g)
   if(tempa_g%l_stk.eq.adss_o) temp_g%l_stk=adss_n
   if(tempa_g%r_stk.eq.adss_o) temp_g%r_stk=adss_n
  end if
 end if

! The following operation updates the grid-map.
 i0=temp%g_cell%x_idx; j0=temp%g_cell%y_idx

!  if(i0.eq.1.and.j0.eq.1) then
!   print*, tempa%address%cv_nmb, tempa%address%idx
!   print*, temp%address%cv_nmb, temp%address%idx
!  end if

 do i=1,4
  if(updated(i0,j0,i).eq.'no ') then
   if(ggd_cell(i0,j0)%ccpt(i)%address.eq.adss_o) then
    ggd_cell(i0,j0)%ccpt(i)%address=adss_n
	updated(i0,j0,i)='yes'
   end if
  end if
 end do
! Variable 'ggd_cell' should be updated due to the possible 
! changes of addresses in the transform from 'tempa' to 'temp'.

! call check_list(cvv,'down  ')
 
 current=>temp

 tempa=>tempa%next
 if(associated(tempa)) then
  head_switch=(tempa%address%idx.ne.head_mark.and.(tempa%address%idx.ne.end_mark))
 else
  exit
 end if
end do

if(cvv(cv_nmb)%cv_type.eq.'circular') then
 cvv(cv_nmb)%begin%next%previous=>cvv(cv_nmb)%eend%previous
 cvv(cv_nmb)%eend%previous%next=>cvv(cv_nmb)%begin%next
end if

deallocate(updated)
! call check_list(cvv,'down  ')
end subroutine form_dc

end module form_discurves