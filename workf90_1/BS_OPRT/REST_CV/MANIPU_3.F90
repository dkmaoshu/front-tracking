module manipulation_3
! This module is for fixing one-side peak occured between two 
! neighboring auxiliary $xy$-type discontinuity cells. The other 
! neighboring diconsintuity cells of both the $xy$-type discontinuity
! cells must be of either $xx$- or $yy$-type.

use discontinuity_curves

implicit none

public  manipulate_3

private fix_peaks, fix_peaks_1


contains


subroutine manipulate_3 

implicit none

if(cvv(1)%status.eq.'asleep ') then
 print*, 'The curve list does not exist.'; return
end if

call fix_peaks

call fix_peaks_1

end subroutine manipulate_3


subroutine fix_peaks

implicit none

type(critical_cell), pointer :: temp, temp_n
type(geo_info) :: g_cell, g_cell_n
type(phy_info) :: p_cell, p_cell_n
integer :: head_mark
logical :: head_switch

logical :: peak, peak_n
real(8) :: distance, distance_n, ratio, ratio_n
type(state) :: difference
integer :: single
character*6 :: side
  
if(associated(cvv(1)%begin%next)) then
 temp=>cvv(1)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if
  
do while(associated(temp).and.head_switch) 

 g_cell=temp%g_cell; p_cell=temp%p_cell
 temp_n=>temp%next
 g_cell_n=temp_n%g_cell; p_cell_n=temp_n%p_cell
  
 if(g_cell%g_type.eq.'xy '.and.g_cell_n%g_type.eq.'xy ') then

  if(temp%previous%g_cell%g_type.eq.'xy ') goto 10 
  if(temp_n%next%g_cell%g_type.eq.'xy ') goto 10 
  if(g_cell%edge(1).ne.g_cell_n%edge(2)) goto 10

  call if_peak(g_cell,1,peak,distance)
  call if_peak(g_cell_n,2,peak_n,distance_n)
   
  if(peak.or.peak_n) then

   call find_single(g_cell,single)
   side=g_cell%point(single)
   call find_single(g_cell_n,single)
   if(g_cell_n%point(single).ne.side) call error_message
   
   if(side.eq.'left  ') then
    difference=temp%p_cell%or_state-temp%p_cell%r_state
    difference=difference+temp_n%p_cell%or_state-temp_n%p_cell%r_state
    if(distance.lt.1.0d-2.and.distance_n.lt.1.0d-2) then
!     ratio=-0.5d0; ratio_n=-0.5d0
     temp%p_cell%or_state=temp%p_cell%r_state+(-1.0d-1)*difference
     temp_n%p_cell%or_state=temp_n%p_cell%r_state+(-1.0d-1)*difference
     temp%previous%p_cell%or_state=temp%previous%p_cell%or_state &
     +0.6d0*difference
     temp_n%next%p_cell%or_state=temp_n%next%p_cell%or_state &
     +0.6d0*difference
    else
!     ratio=distance_n/(distance+distance_n)
!     ratio_n=distance/(distance+distance_n)
     temp%p_cell%or_state=temp%p_cell%r_state+0.5d0*difference
     temp_n%p_cell%or_state=temp_n%p_cell%r_state+0.5d0*difference
    end if
   else
    difference=temp%p_cell%or_state-temp%p_cell%l_state
    difference=difference+temp_n%p_cell%or_state-temp_n%p_cell%l_state
    if(distance.lt.1.0d-2.and.distance_n.lt.1.0d-2) then
     temp%p_cell%or_state=temp%p_cell%l_state+(-1.0d-1)*difference
     temp_n%p_cell%or_state=temp_n%p_cell%l_state+(-1.0d-1)*difference
     temp%previous%p_cell%or_state=temp%previous%p_cell%or_state &
     +0.6d0*difference
     temp_n%next%p_cell%or_state=temp_n%next%p_cell%or_state &
     +0.6d0*difference
    else
     temp%p_cell%or_state=temp%p_cell%l_state+0.5d0*difference
     temp_n%p_cell%or_state=temp_n%p_cell%l_state+0.5d0*difference
    end if
   end if
  end if   
 end if

 10 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 endif
  
end do   


contains


subroutine if_peak(g_cell,num,peak,distance)

implicit none

type(geo_info), intent(in) :: g_cell
integer, intent(in) :: num
logical, intent(out) :: peak
real(8), intent(out) :: distance

integer :: single
real(8) :: other_distance

peak=.false.
call find_single(g_cell,single)
  
select case(single)
  
 case(1) 
  distance=dabs(g_cell%dis_posi(num)+0.5d0)
  other_distance=dabs(g_cell%dis_posi(3-num)+0.5d0)
 
 case(2)
  if(g_cell%edge(num).eq.1) then
   distance=dabs(0.5d0-g_cell%dis_posi(num))
   other_distance=dabs(g_cell%dis_posi(3-num)+0.5d0)
  else
   distance=dabs(g_cell%dis_posi(num)+0.5d0)
   other_distance=dabs(0.5d0-g_cell%dis_posi(3-num)) 
  end if  
  
 case(3) 
  distance=dabs(0.5-g_cell%dis_posi(num))
  other_distance=dabs(0.5-g_cell%dis_posi(3-num))
 
 case(4)
  if(g_cell%edge(num).eq.4) then
   distance=dabs(0.5d0-g_cell%dis_posi(num))
   other_distance=dabs(g_cell%dis_posi(3-num)+0.5d0)
  else
   distance=dabs(g_cell%dis_posi(num)+0.5d0)
   other_distance=dabs(0.5d0-g_cell%dis_posi(3-num)) 
  end if          

end select

if(distance.lt.1.0d-2) peak=.true.

if(distance.lt.1.0d-1.and.other_distance.gt.1.0d-1) then
 if(distance/other_distance.lt.6.0d0) peak=.true. 
end if

end subroutine if_peak


end subroutine fix_peaks


subroutine fix_peaks_1

implicit none

type(critical_cell), pointer :: temp, temp_p, temp_n
type(geo_info) :: g_cell, g_cell_p, g_cell_n
type(phy_info) :: p_cell, p_cell_p, p_cell_n
integer :: head_mark
logical :: head_switch

logical :: peak, peak_p, peak_n
type(state) :: difference
integer::single
character*6 :: side, other_side
  
if(associated(cvv(1)%begin%next)) then
 temp=>cvv(1)%begin%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if
  
do while(associated(temp).and.head_switch) 

 temp_p=>temp%previous; temp_n=>temp%next
 g_cell=temp%g_cell; g_cell_p=temp_p%g_cell; g_cell_n=temp_n%g_cell
  
 if(g_cell%g_type.eq.'xy '.and.g_cell_p%g_type.eq.'xy ') then
  if(g_cell_n%g_type.eq.'xy ') then

   if(g_cell%edge(1).ne.g_cell_n%edge(2)) goto 10
   if(g_cell_p%edge(1).ne.g_cell%edge(2)) goto 10   

   call if_peak(g_cell,peak)
   call if_peak(g_cell_p,peak_p)
   call if_peak(g_cell_n,peak_n)
   
   if(peak.or.peak_p.or.peak_n) then
    
	call find_single(g_cell,single)
    side=g_cell%point(single)

    call find_single(g_cell_p,single)
	if(g_cell_p%point(single).ne.side) call error_message
    call find_single(g_cell_n,single)
	if(g_cell_n%point(single).ne.side) call error_message

    if(side.eq.'left  ') then
     difference=temp%p_cell%or_state-temp%p_cell%r_state
	 difference=difference+temp_p%p_cell%or_state-temp_p%p_cell%r_state
     difference=difference+temp_n%p_cell%or_state-temp_n%p_cell%r_state
     temp%p_cell%or_state=temp%p_cell%r_state+(-1.0d-1)*difference
     temp_p%p_cell%or_state=temp_p%p_cell%r_state+(-1.0d-1)*difference
     temp_n%p_cell%or_state=temp_n%p_cell%r_state+(-1.0d-1)*difference
    else
     difference=temp%p_cell%or_state-temp%p_cell%l_state
	 difference=difference+temp_p%p_cell%or_state-temp_p%p_cell%l_state
     difference=difference+temp_n%p_cell%or_state-temp_n%p_cell%l_state
     temp%p_cell%or_state=temp%p_cell%l_state+(-1.0d-1)*difference
     temp_p%p_cell%or_state=temp_p%p_cell%l_state+(-1.0d-1)*difference
     temp_n%p_cell%or_state=temp_n%p_cell%l_state+(-1.0d-1)*difference
	end if

    temp_p%previous%p_cell%or_state=temp_p%previous%p_cell%or_state &
    +0.65d0*difference
    temp_n%next%p_cell%or_state=temp_n%next%p_cell%or_state &
    +0.65d0*difference
   
   end if
  end if   
 end if

 10 temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 endif
  
end do   


contains


subroutine if_peak(g_cell,peak)

implicit none

type(geo_info), intent(in) :: g_cell
logical, intent(out) :: peak

real(8) :: distance
integer :: single
real(8) :: other_distance

peak=.false.
call find_single(g_cell,single)
  
select case(single)
  
 case(1) 
  distance=dabs(g_cell%dis_posi(1)+0.5d0)
  other_distance=dabs(g_cell%dis_posi(2)+0.5d0)
 
 case(2)
  if(g_cell%edge(1).eq.1) then
   distance=dabs(0.5d0-g_cell%dis_posi(1))
   other_distance=dabs(g_cell%dis_posi(2)+0.5d0)
  else
   distance=dabs(g_cell%dis_posi(1)+0.5d0)
   other_distance=dabs(0.5d0-g_cell%dis_posi(2)) 
  end if  
  
 case(3) 
  distance=dabs(0.5-g_cell%dis_posi(1))
  other_distance=dabs(0.5-g_cell%dis_posi(2))
 
 case(4)
  if(g_cell%edge(1).eq.4) then
   distance=dabs(0.5d0-g_cell%dis_posi(1))
   other_distance=dabs(g_cell%dis_posi(2)+0.5d0)
  else
   distance=dabs(g_cell%dis_posi(1)+0.5d0)
   other_distance=dabs(0.5d0-g_cell%dis_posi(2)) 
  end if          

end select

if(distance.lt.1.0d-2.or.other_distance.lt.1.0d-2) peak=.true.

end subroutine if_peak


end subroutine fix_peaks_1


end module manipulation_3


