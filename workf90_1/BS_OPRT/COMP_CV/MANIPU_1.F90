module manipulation_1
! This module is for the following manipulation in the computation. When
! there are two neighboring single $xy$-type discontinuity cells on the
! auxiliary discontinuity curve and their ordinary cell-averages (maybe 
! only one of them) exceed the left and right bounds, the exceeded parts
! of the ordinary cell-averages should go to their neighboring $xx$- or
! $yy$-type discontinuity cells. Also when there are sigle $xy$-type
! discontinuity cells and the two neighboring $xx$- or $yy$-type 
! discontiuity cells are both sigle, the ordinary cell-averages of
! the single $xy$-type discontinuity cells are also modified and the 
! exceeded amount of the physical quantity goes evenly to the two
! neighbors.

use auxiliary_discontinuity_curves

implicit none

public  manipulate_1

private truncation_1, distribution_1 


contains


subroutine manipulate_1 

implicit none

integer :: ii

do ii=1, cvn
 if(acvv(1)%status.eq.'awake  ') then
  call truncation_1(ii)
  call distribution_1(ii)
 end if 
end do

end subroutine manipulate_1


subroutine truncation_1(cv_nmb)

implicit none

integer, intent(in) :: cv_nmb

type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell, g_cell_p, g_cell_n
type(phy_info) :: p_cell
type(state) :: difference
logical :: if_truncate
integer :: head_mark
logical :: head_switch
  
if(associated(acvv(cv_nmb)%begin%next)) then
 temp=>acvv(cv_nmb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if
  
do while(associated(temp).and.head_switch) 
  
 call clean_up_g_cell(g_cell); call clean_up_g_cell(g_cell_n)
 call clean_up_g_cell(g_cell_p)
 call clean_up_p_cell(p_cell) 
  
 g_cell=temp%g_cell
 g_cell_n=temp%next%g_cell; g_cell_p=temp%previous%g_cell
 if_truncate=.false.
 temp%memo=0.0d0
  
 if(g_cell%g_type.eq.'xy '.and.g_cell_n%g_type.eq.'xy ') then 
  if_truncate=.true.
 end if

 if(g_cell%g_type.eq.'xy '.and.g_cell_p%g_type.eq.'xy ') then
  if_truncate=.true.
 end if
  
 if(g_cell%g_type.eq.'xy ') then
  if(g_cell_p%g_type.ne.'xy '.and.temp%previous%partner.eq.'single') then
   if(g_cell_n%g_type.ne.'xy '.and.temp%next%partner.eq.'single') then  
    if_truncate=.true.
   end if
  end if  
 end if
  
 if(g_cell%g_type.eq.'xx '.or.g_cell%g_type.eq.'yy ') then
  if(temp%partner.eq.'single') then
   if(g_cell_p%g_type.eq.'xy '.or.g_cell_n%g_type.eq.'xy ') then
    if_truncate=.true.
   end if
  end if
 end if
  
 if(if_truncate) then
  p_cell=temp%p_cell
  call truncation(g_cell,p_cell,difference)
  temp%p_cell=p_cell
  temp%memo=difference
  
!   if(dabs(difference%value).gt.1.0d-6) then
!    write(*,*) 'TRUNCATED!!!', g_cell%x_idx, g_cell%y_idx 
!   end if

 end if   
  
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 endif
  
end do   


contains


subroutine truncation(g_cell,p_cell,difference)

implicit none

type(geo_info), intent(in) :: g_cell
type(phy_info), intent(inout) :: p_cell
type(state), intent(out) :: difference

integer :: single
type(state) :: uuu

difference=0.0d0
select case(g_cell%g_type)
  
 case('xx ')
    
  if(p_cell%l_state.gt.p_cell%r_state) then
   if(p_cell%or_state.gt.p_cell%l_state) then
    uuu=(1.0d0-1.0d-5)*p_cell%l_state+1.0d-5*p_cell%r_state
    difference=p_cell%or_state-uuu
    p_cell%or_state=uuu
   else
    if(p_cell%or_state.lt.p_cell%r_state) then
     uuu=(1.0d0-1.0d-5)*p_cell%r_state+1.0d-5*p_cell%l_state
	 difference=p_cell%or_state-uuu
     p_cell%or_state=uuu
    end if
   end if	     
  else
   if(p_cell%or_state.gt.p_cell%r_state) then
    uuu=(1.0d0-1.0d-5)*p_cell%r_state+1.0d-5*p_cell%l_state
    difference=p_cell%or_state-uuu
    p_cell%or_state=uuu
   else
    if(p_cell%or_state.lt.p_cell%l_state) then
     uuu=(1.0d0-1.0d-5)*p_cell%l_state+1.0d-5*p_cell%r_state
     difference=p_cell%or_state-uuu
     p_cell%or_state=uuu
    end if 
   end if
  end if  
  
 case('yy ')
  
  if(p_cell%l_state.gt.p_cell%r_state) then
   if(p_cell%or_state.gt.p_cell%l_state) then
    uuu=(1.0d0-1.0d-5)*p_cell%l_state+1.0d-5*p_cell%r_state
    difference=p_cell%or_state-uuu
    p_cell%or_state=uuu
   else
    if(p_cell%or_state.lt.p_cell%r_state) then
     uuu=(1.0d0-1.0d-5)*p_cell%r_state+1.0d-5*p_cell%l_state
	 difference=p_cell%or_state-uuu
     p_cell%or_state=uuu
    end if
   end if	     
  else
   if(p_cell%or_state.gt.p_cell%r_state) then
    uuu=(1.0d0-1.0d-5)*p_cell%r_state+1.0d-5*p_cell%l_state
    difference=p_cell%or_state-uuu
    p_cell%or_state=uuu
   else
    if(p_cell%or_state.lt.p_cell%l_state) then
     uuu=(1.0d0-1.0d-5)*p_cell%l_state+1.0d-5*p_cell%r_state
     difference=p_cell%or_state-uuu
     p_cell%or_state=uuu
    end if 
   end if
  end if  
   
 case('xy ')

  call find_single(g_cell,single)
  
  select case(g_cell%point(single))
    
   case('left  ')
   
    if(p_cell%l_state.gt.p_cell%r_state) then
     if(p_cell%or_state.gt.0.5d0*(p_cell%l_state+p_cell%r_state)) then
      difference=p_cell%or_state-0.5d0*(p_cell%l_state+p_cell%r_state)
      p_cell%or_state=0.5d0*(p_cell%l_state+p_cell%r_state)
     else
      if(p_cell%or_state.lt.p_cell%r_state) then
       uuu=(1.0d0-1.0d-5)*p_cell%r_state+1.0d-5*p_cell%l_state
       difference=p_cell%or_state-uuu
       p_cell%or_state=uuu
      end if
     end if	     
    else
     if(p_cell%or_state.gt.p_cell%r_state) then
      uuu=(1.0d0-1.0d-5)*p_cell%r_state+1.0d-5*p_cell%l_state
      difference=p_cell%or_state-uuu
      p_cell%or_state=uuu
     else
      if(p_cell%or_state.lt.0.5d0*(p_cell%l_state+p_cell%r_state)) then
       difference=p_cell%or_state-0.5d0*(p_cell%l_state+p_cell%r_state)
       p_cell%or_state=0.5d0*(p_cell%l_state+p_cell%r_state)
      end if 
     end if
    end if  
    
   case('right ')
   
    if(p_cell%l_state.gt.p_cell%r_state) then
     if(p_cell%or_state.gt.p_cell%l_state) then
      uuu=(1.0d0-1.0d-5)*p_cell%l_state+1.0d-5*p_cell%r_state
      difference=p_cell%or_state-uuu
      p_cell%or_state=uuu
     else
      if(p_cell%or_state.lt.0.5d0*(p_cell%l_state+p_cell%r_state)) then
       difference=p_cell%or_state-0.5d0*(p_cell%l_state+p_cell%r_state)
       p_cell%or_state=0.5d0*(p_cell%l_state+p_cell%r_state)
      end if
     end if	     
    else
     if(p_cell%or_state.gt.0.5d0*(p_cell%l_state+p_cell%r_state)) then
      difference=p_cell%or_state-0.5d0*(p_cell%l_state+p_cell%r_state)
      p_cell%or_state=0.5d0*(p_cell%l_state+p_cell%r_state)
     else
      if(p_cell%or_state.lt.p_cell%l_state) then
       uuu=(1.0d0-1.0d-5)*p_cell%l_state+1.0d-5*p_cell%r_state
       difference=p_cell%or_state-uuu
       p_cell%or_state=uuu
      end if 
     end if
    end if  
     
   case default; call error_message

  end select

end select


end subroutine truncation


end subroutine truncation_1
   

subroutine distribution_1(cv_nmb)

implicit none

integer, intent(in) :: cv_nmb

type(auxi_crit_cell), pointer :: temp
type(geo_info) :: g_cell, g_cell_n, g_cell_p
type(phy_info) :: p_cell, p_cell_p, p_cell_n
type(state) :: difference
logical :: if_truncate
integer :: head_mark
logical :: head_switch
  
if(associated(acvv(cv_nmb)%begin%next)) then
 temp=>acvv(cv_nmb)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if
  
do while(associated(temp).and.head_switch) 

 call clean_up_g_cell(g_cell); call clean_up_g_cell(g_cell_p)
 call clean_up_g_cell(g_cell_n)
 call clean_up_p_cell(p_cell); call clean_up_p_cell(p_cell_n)
 call clean_up_p_cell(p_cell_p)

 g_cell=temp%g_cell
 g_cell_p=temp%previous%g_cell; g_cell_n=temp%next%g_cell
 p_cell_n=temp%next%p_cell; p_cell_p=temp%previous%p_cell
 difference=temp%memo
 
! if(g_cell%g_type.eq.'xy '.and.g_cell_n%g_type.eq.'xy ') then
!  p_cell_p%or_state=p_cell_p%or_state+difference
!  if(temp%previous%partner.eq.'previous') then
!   temp%previous%previous%p_cell%or_state= &
!   temp%previous%previous%p_cell%or_state+difference
!  end if
! else
!  if(g_cell%g_type.eq.'xy '.and.g_cell_p%g_type.eq.'xy ') then
!   p_cell_n%or_state=p_cell_n%or_state+difference
!   if(temp%next%partner.eq.'next    ') then
!    temp%next%next%p_cell%or_state= &
!    temp%next%next%p_cell%or_state+difference
!   end if
!  else
   p_cell_p%or_state=p_cell_p%or_state+0.5d0*difference
   if(temp%previous%partner.eq.'previous') then
    temp%previous%previous%p_cell%or_state= &
    temp%previous%previous%p_cell%or_state+0.5d0*difference
   end if
   p_cell_n%or_state=p_cell_n%or_state+0.5d0*difference
   if(temp%next%partner.eq.'next    ') then
    temp%next%next%p_cell%or_state= &
    temp%next%next%p_cell%or_state+0.5d0*difference
   end if     
!  end if
! end if

 temp%previous%p_cell=p_cell_p
 temp%next%p_cell=p_cell_n
 
 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 endif
  
end do   

end subroutine distribution_1


end module manipulation_1