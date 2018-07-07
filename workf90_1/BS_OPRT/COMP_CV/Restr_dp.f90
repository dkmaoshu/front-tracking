module restriction_of_positions

use	auxiliary_discontinuity_curves
! 'auxi_cv.f90'

type(auxi_crit_cell), pointer :: temp, boundary_end
type(geo_info) :: g_cell

public restrict_positions
private	temp, boundary_end, g_cell


contains


subroutine restrict_positions(step)

implicit none
integer, intent(in) :: step

integer :: i

do i=1, cvn
 if(acvv(i)%status.eq.'awake') then
  call restrict_posi(acvv(i))
 end if
end do


contains


subroutine restrict_posi(acv)

implicit none

type(auxi_discv), intent(inout) :: acv

type(adss_info) :: address
integer :: head_mark, end_mark
logical :: head_switch
integer :: i, difference
character*8 :: partner
type(geo_info) :: gn_cell, gp, gpn, gn, gnn
real(8) :: restriction_1, restriction_2

select case(step)
 case(0)
  restriction_1=0.1d0; restriction_2=0.1d0
 case(1)
  restriction_1=0.5d0; restriction_2=0.45d0
 case default; call error_message
end select

if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
 boundary_end=>acv%eend_boundary%previous
 if(acv%eend%previous%address.eq.boundary_end%address) then
  end_mark=error_index
 else
  end_mark=acv%eend%previous%next%address%idx
 end if

 do while(associated(temp).and.head_switch) 
  difference=error_index
  address=temp%address; partner=temp%partner
  g_cell=temp%g_cell

! Restrict discontinuity positions in the current critical cell.
  if(partner.eq.'single  ') then
! The case of no-partner.   
   do i=1,2
    if(g_cell%dis_posi(i).gt.0.d0) then
     if(g_cell%dis_posi(i).gt.0.5d0+restriction_1) then
      g_cell%dis_posi(i)=0.5d0+restriction_1-1.0d-6
      print*, 'Restricted', g_cell%x_idx, g_cell%y_idx
     end if
    else 
     if(g_cell%dis_posi(i).lt.-0.5d0-restriction_1) then
      g_cell%dis_posi(i)=-0.5d0-restriction_1+1.0d-6
      print*, 'Restricted', g_cell%x_idx, g_cell%y_idx
     end if
    end if										
   end do
  else
! the case of with-partner.
   if(partner.eq.'next  ') then
    gn_cell=temp%next%g_cell
    difference=g_cell%x_idx-gn_cell%x_idx+g_cell%y_idx-gn_cell%y_idx
    do i=1,2
     if(difference.gt.0) then
      if(g_cell%dis_posi(i).lt.-1.5d0-restriction_2) then
       g_cell%dis_posi(i)=-1.5d0-restriction_2+1.0d-6
       print*, 'Restricted', g_cell%x_idx, g_cell%y_idx
      endif
      if(g_cell%dis_posi(i).gt.0.5d0+restriction_2) then
       g_cell%dis_posi(i)=0.5d0+restriction_2+1.0d-6
       print*, 'Restricted', g_cell%x_idx, g_cell%y_idx
      end if
     end if
     if(difference.lt.0) then
      if(g_cell%dis_posi(i).lt.-0.5d0-restriction_2) then
       g_cell%dis_posi(i)=-0.5d0-restriction_2+1.0d-6
       print*, 'Restricted', g_cell%x_idx, g_cell%y_idx
      end if
      if(g_cell%dis_posi(i).gt.1.5d0+restriction_2) then
       g_cell%dis_posi(i)=1.5d0+restriction_2-1.0d-6
       print*, 'Restricted', g_cell%x_idx, g_cell%y_idx
      end if
     end if
    end do
    do i=1,2
     gn_cell%dis_posi(i)=g_cell%dis_posi(i)+difference
    end do
    temp%next%g_cell=gn_cell
   end if
  end if

! Update the discontinuity position in the previous next-door 
! critical cell.
  if(associated(temp%p_nxt_dr)) then
   gp=temp%p_nxt_dr%g_cell
   if(g_cell%edge(1).eq.1.or.g_cell%edge(1).eq.3) then
    difference=g_cell%x_idx-gp%x_idx
   else
    difference=g_cell%y_idx-gp%y_idx
   end if 
   gp%dis_posi(2)=g_cell%dis_posi(1)+difference
   if(temp%p_nxt_dr%partner.ne.'single  ') then
    if(temp%p_nxt_dr%partner.eq.'previous') then
     gpn=temp%p_nxt_dr%previous%g_cell
    else
     gpn=temp%p_nxt_dr%next%g_cell
    end if
    difference=gp%x_idx-gpn%x_idx+gp%y_idx-gpn%y_idx
    gpn%dis_posi(2)=gp%dis_posi(2)+difference
    if(temp%p_nxt_dr%partner.eq.'previous') then
     temp%p_nxt_dr%previous%g_cell=gpn
    else
     temp%p_nxt_dr%next%g_cell=gpn
    end if
   end if
   temp%p_nxt_dr%g_cell=gp
  end if 

! Update the discontinuity position in the next next-door
! critical cell.
  if(associated(temp%n_nxt_dr)) then
   gn=temp%n_nxt_dr%g_cell
   if(g_cell%edge(2).eq.1.or.g_cell%edge(2).eq.3) then
    difference=g_cell%x_idx-gn%x_idx
   else
    difference=g_cell%y_idx-gn%y_idx
   end if
   gn%dis_posi(1)=g_cell%dis_posi(2)+difference
   if(temp%n_nxt_dr%partner.ne.'single  ') then
    if(temp%n_nxt_dr%partner.eq.'previous') then
     gnn=temp%n_nxt_dr%previous%g_cell
    else
     gnn=temp%n_nxt_dr%next%g_cell
    end if
    difference=gn%x_idx-gnn%x_idx+gn%y_idx-gnn%y_idx
    gnn%dis_posi(1)=gn%dis_posi(1)+difference
    if(temp%n_nxt_dr%partner.eq.'previous') then
     temp%n_nxt_dr%previous%g_cell=gnn
    else
     temp%n_nxt_dr%next%g_cell=gnn
    end if
   end if
   temp%n_nxt_dr%g_cell=gn
  end if
   
  temp%g_cell=g_cell

! Go to the next critical cell on the curve.
  temp=>temp%next
  if(associated(temp)) then
   head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
  else
   exit
  end if
 end do
end if

end subroutine restrict_posi


end subroutine restrict_positions

end module restriction_of_positions