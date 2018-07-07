module auxiliary_curves_production
! This module is for production of auxilary discontinuity curves.

use discontinuity_curves
! 'discontinuity_curve.f90'

use	auxiliary_discontinuity_curves
! 'auxi_cv.f90'

use auxi_cc_prod
! 'auxi_ccp.f90'

use auxiliary_list_fix
! 'isl_fix.f90'

!use recover_mid_positions
!! 'rcv_mp.f90'

 use node_cells
!! 'nd_cls.f90'

! use scan_conservation_on_list
! 'sc_conls.f90'

use clean_up
! 'clean_up.f90'

implicit none

public  production
private produce_auxi_cv, set_nextdoor, produce, list_fixing


contains


subroutine production(acvg_type)

implicit none

character*3, intent(in) :: acvg_type

integer :: i

! type(state) :: fff

allocate(cvv(cvn)); allocate(acvv(cvn)); allocate(ndd(ndn))
call give_cv; call give_nd

! call clean_up_discontinuities

! call check_list_c(1,'down  ')
! call check_ring('counter ')

call change_head_4_circular

! acvg_type='yy '

do i=1,cvn
 acvv(i)%status=cvv(i)%status
 if(cvv(i)%status.eq.'awake '.or.cvv(i)%status.eq.'yawn  ') then
  acvv(i)%cv_type=cvv(i)%cv_type
  acvv(i)%begin_end=cvv(i)%begin_end
  acvv(i)%end_end=cvv(i)%end_end
  acvv(i)%total=cvv(i)%total; acvv(i)%wave=cvv(i)%wave
  call produce_auxi_cv(acvg_type,i)
!   call check_list_ac(1,'down  ')
  call set_nextdoor(acvv(i))

 end if
end do

! fff=acvv(1)%begin%c_cell%l_state(0)
! call check_list_ac(1,'down  ')

call list_fixing

! call check_list_ac(1,'down  ')

!call smoothen_xy

! call check_partner_conserv_ac(1)

do i=1,cvn
 call deletee(cvv(i))
end do

call curves_type_change('auxila')

! call delete(acvv(8))

call get_cv
call get_acv
call get_nd

deallocate(cvv); deallocate(acvv); deallocate(ndd)

end subroutine production


subroutine produce_auxi_cv(acvg_type,cv_nmb)

character*3, intent(in) :: acvg_type
integer, intent(in) :: cv_nmb

type(critical_cell), pointer :: temp
type(auxi_crit_cell), pointer :: tempa, currenta
type(geo_info) :: ag_c(2), g_c, g_a, g_aa
type(phy_info) :: ap_c(2), p_c, p_a
integer :: k, steps, head_mark
logical :: head_switch

call creat_acv(cv_nmb)
if(cvv(cv_nmb)%status.eq.'yawn  ') return

currenta=>acvv(cv_nmb)%begin
if(associated(cvv(cv_nmb)%begin%next)) then
 temp=>cvv(cv_nmb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.

 g_c=temp%g_cell; call clean_up_g_cell(g_a); call clean_up_g_cell(g_aa)
 if(cvv(cv_nmb)%cv_type.eq.'circular') then
! If the discontinuity curve is circular, production should start
! with an appropriate critical cell, either with a one of either 
! 'xx' or 'yy' type, or with a one of 'xy' type that will have 
! partnership with its next critical cell.
  if(temp%g_cell%g_type.eq.'xy') then
   g_a=temp%previous%g_cell
   call with_geo(g_c,g_a,ag_c(1),'previous')
   if(ag_c(1)%g_type.eq.'xx'.or. ag_c(1)%g_type.eq.'yy') then
    g_aa=temp%previous%previous%g_cell
    call with_geo(g_a,g_aa,ag_c(2),'previous')
    if(ag_c(2)%g_type.ne.acvg_type) then
     temp=>temp%previous; head_mark=temp%address%idx 
    end if
   end if
  end if
 end if      

! To produce the auxiliary discontinuity curve.
 do while(associated(temp).and.head_switch) 
! To prepare the necessary data.
  call clean_up_g_cell(g_a); call clean_up_g_cell(g_aa); call clean_up_p_cell(p_a)
  g_c=temp%g_cell; p_c=temp%p_cell
  if(associated(temp%next)) then
   g_a=temp%next%g_cell
   p_a=temp%next%p_cell
  end if
! To decide the case (steps).
  if(g_c%g_type.ne.'xy') then
   steps=1
  else
   if(associated(temp%next)) then
    call with_geo(g_c,g_a,ag_c(1),'next    ')
    if(ag_c(1)%g_type.eq.acvg_type) then
     steps=2
    else
     if(ag_c(1)%g_type.eq.'xx'.or.ag_c(1)%g_type.eq.'yy') then
      if(associated(temp%next%next)) then
       g_aa=temp%next%next%g_cell
       call with_geo(g_a,g_aa,ag_c(2),'next    ')
       if(ag_c(2)%g_type.ne.acvg_type) then
        steps=2
       else
        steps=1
       endif
      else
       steps=2
      end if
     else
      steps=1
     end if
    end if
   else
    steps=1
   end if
  end if
! To produce the geometry and physics cells for production.
  call produce(g_c,g_a,p_c,p_a,ag_c,ap_c,steps)
! To creat the auxiliary critical cell and link it to the curve. 
  do k=1, steps
   allocate(tempa)
   nullify(tempa%previous); nullify(tempa%next)
   call clean_up_c_cell(tempa%c_cell)

   tempa%g_cell=ag_c(k); tempa%p_cell=ap_c(k)
   tempa%address=temp%address
   tempa%temp=>temp
   tempa%l_stk=temp%l_stk; tempa%r_stk=temp%r_stk
   tempa%c_cell%l_state(0)=ap_c(k)%l_state
   tempa%c_cell%r_state(0)=ap_c(k)%r_state
   tempa%c_cell%left_room_for_swap=ap_c(k)%l_state
   tempa%c_cell%right_room_for_swap=ap_c(k)%r_state
   tempa%c_cell%or_state(0)=ap_c(k)%or_state
   tempa%updated='no '

   if(steps.eq.1) then
    tempa%partner='single  '
   else
    if(k.eq.1) then
     tempa%partner='next    '
    else
     tempa%partner='previous'
    end if
   end if  		 		 	    

!   call insert(currenta,tempa,'after ')
   currenta%next=>tempa
   if(temp%address%idx.ne.head_mark) tempa%previous=>currenta
   currenta=>tempa
   
   temp=>temp%next
   if(associated(temp)) then
    head_switch=(temp%address%idx.ne.head_mark)
   else
    exit
   end if

  end do

 end do

 acvv(cv_nmb)%eend%previous=>currenta 
 if(acvv(cv_nmb)%cv_type.eq.'circular') then
  currenta%next=>acvv(cv_nmb)%begin%next
  acvv(cv_nmb)%begin%next%previous=>currenta
 end if

end if

acvv(cv_nmb)%begin_boundary%next=>acvv(cv_nmb)%begin%next
acvv(cv_nmb)%eend_boundary%previous=>acvv(cv_nmb)%eend%previous

! call visit(adss_info(1,32),temp)
  
end subroutine produce_auxi_cv


subroutine set_nextdoor(acvv)
! Set the next-doorship for the auiliary critical cells on the 
! curve.

implicit none

type(auxi_discv) :: acvv
type(auxi_crit_cell), pointer :: tempa
integer :: head_mark
logical :: head_switch

if(acvv%status.eq.'yawn  ') return
tempa=>acvv%begin%next
head_mark=tempa%address%idx; head_switch=.true.
do while(associated(tempa).and.head_switch) 
 nullify(tempa%p_nxt_dr); nullify(tempa%n_nxt_dr)

 if(associated(tempa%previous)) then

  if(tempa%partner.eq.'single  ') then
   tempa%p_nxt_dr=>tempa%previous
   tempa%previous%n_nxt_dr=>tempa
   if(tempa%previous%partner.eq.'previous') then
    tempa%previous%previous%n_nxt_dr=>tempa
   endif
  else
   if(tempa%partner.eq.'previous') then
    tempa%p_nxt_dr=>tempa%previous%previous
   else
    tempa%p_nxt_dr=>tempa%previous
    tempa%previous%n_nxt_dr=>tempa
    if(tempa%previous%partner.eq.'previous') then
     tempa%previous%previous%n_nxt_dr=>tempa
    end if      	 	  	   
   end if
  end if	 
 
 end if

 tempa=>tempa%next
 if(associated(tempa)) then
  head_switch=(tempa%address%idx.ne.head_mark)
 else
  exit
 end if

 if(acvv%cv_type.eq.'circular') then
  if(tempa%partner.eq.'single') then
   tempa%p_nxt_dr=>tempa%previous
   tempa%previous%n_nxt_dr=>tempa
   if(tempa%previous%partner.eq.'previous') then
    tempa%previous%previous%n_nxt_dr=>tempa
   endif
  else
   if(tempa%partner.eq.'previous') then
    tempa%p_nxt_dr=>tempa%previous%previous
   else
    tempa%p_nxt_dr=>tempa%previous
    tempa%previous%n_nxt_dr=>tempa
    if(tempa%previous%partner.eq.'previous') then
     tempa%previous%previous%n_nxt_dr=>tempa
    end if      	 	  	   
   end if
  end if
 end if

end do

end subroutine set_nextdoor


subroutine produce(g_c,g_a,p_c,p_a,ag_c,ap_c,steps)
! To proiduce the geometry and physics cells.
implicit none
type(geo_info) :: ag_c(2), g_c, g_a 
type(phy_info) :: ap_c(2), p_c, p_a
integer :: steps
intent(in) :: g_c, g_a, p_c, p_a, steps
intent(out) :: ag_c, ap_c

if(steps.eq.1) then
 ag_c(1)=g_c; ap_c(1)=p_c
else
 call with_geo(g_c,g_a,ag_c(1),'next    ')
 call with_phy(g_c,g_a,p_c,p_a,ap_c(1))
 call with_geo(g_a,g_c,ag_c(2),'previous')
 call with_phy(g_a,g_c,p_a,p_c,ap_c(2))
end if

end subroutine produce

end module auxiliary_curves_production