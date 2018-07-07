module manipulation_2
! When two stacked discontinuity cells have their discontinuity positions
! propagate across to each other, this module fixes the crossing. In 
! doing so, we deny any surface tension, either the real or numerical 
! ones.

use auxiliary_discontinuity_curves

implicit none

public  manipulate_2


contains


subroutine manipulate_2 

implicit none

if(acvv(1)%status.eq.'asleep ') then
 print*, 'The curve list does not exist.'
 return
end if

call adjust_crossing(acvv(1))

end subroutine manipulate_2


subroutine adjust_crossing(acv)

implicit none

type(auxi_discv) ::acv
type(auxi_crit_cell), pointer :: temp, tempn_l, tempn_r
type(geo_info) :: g_cell, g_celln_l, g_celln_r
integer :: head_mark
logical :: head_switch

logical :: left_cross, right_cross
integer :: numn_l, numn_r
real(8) :: xx
  
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
else
 print*, 'The curve list is empty.'; return
end if
  
do while(associated(temp).and.head_switch) 

 if(temp%l_stk.eq.adss_info(0,0).and.temp%l_stk.eq.adss_info(0,0)) then
  goto 10
 end if
  
 call clean_up_g_cell(g_cell)
 call clean_up_g_cell(g_celln_l); call clean_up_g_cell(g_celln_r)
  
 g_cell=temp%g_cell

 if(temp%l_stk.ne.adss_info(0,0)) then
  call visit(temp%l_stk,tempn_l)
  g_celln_l=tempn_l%g_cell
  call if_crossing(g_celln_l,'left  ',1,left_cross,numn_l)
 end if
 if(temp%r_stk.ne.adss_info(0,0)) then
  call visit(temp%r_stk,tempn_r)
  g_celln_r=tempn_r%g_cell
  call if_crossing(g_celln_r,'right ',1,right_cross,numn_r)
 end if

 if(left_cross.and.right_cross) call error_message
  
 if(left_cross.or.right_cross) then
  
  if(left_cross) then
   xx=0.5d0*(g_cell%dis_posi(1)+g_celln_l%dis_posi(numn_l))
   g_cell%dis_posi(1)=xx; g_celln_l%dis_posi(numn_l)=xx
   temp%g_cell=g_cell; tempn_l%g_cell=g_celln_l
   call modify_neigh(temp,1) 
   call modify_neigh(tempn_l,numn_l)
  else
   xx=0.5d0*(g_cell%dis_posi(1)+g_celln_r%dis_posi(numn_r))
   g_cell%dis_posi(1)=xx; g_celln_r%dis_posi(numn_r)=xx 
   temp%g_cell=g_cell; tempn_r%g_cell=g_celln_r
   call modify_neigh(temp,1)
   call modify_neigh(tempn_r,numn_r)
  end if
 end if
    
 10  temp=>temp%next
 if(associated(temp)) then
  head_switch=(temp%address%idx.ne.head_mark)
 endif
  
end do   
   
   
contains
   
   
subroutine if_crossing(g_celln,side,num,cross,numn)
   
implicit none
   
type(geo_info), intent(in) :: g_celln
character*6, intent(in) :: side
integer, intent(in) :: num
logical, intent(out) :: cross
integer, intent(out) :: numn

character*6 :: other_side, direction
integer :: single

cross=.false.

! First, determine the possible crossing positions in stacked neighbor.
if(g_celln%edge(1).eq.g_cell%edge(num)) then
 numn=1
else
 if(g_celln%edge(2).eq.g_cell%edge(num)) then
  numn=2
 else
  return
 end if
end if

!Then, determine the direcition of possible crossing.
select case(g_cell%g_type) 

 case('xx ')

  if(g_cell%point(1).eq.side.and.g_cell%point(4).eq.side) then
   direction='left  ' 
  else
   if(g_cell%point(2).eq.side.and.g_cell%point(3).eq.side) then
    direction='right '
   else
    call error_message
   end if
  end if

 case('yy ')

  if(g_cell%point(1).eq.side.and.g_cell%point(2).eq.side) then
   direction='down  ' 
  else
   if(g_cell%point(3).eq.side.and.g_cell%point(4).eq.side) then
    direction='up    '
   else
    call error_message
   end if
  end if

 case('xy ')
  
  call find_single(g_cell,single)        
  
  if(g_cell%point(single).eq.side) then
   
   select case(single)
    
    case(1)
     if(g_cell%edge(num).eq.1) then
      direction='left  '
     else
      if(g_cell%edge(num).eq.4) then
       direction='down  '
      else
       call error_message
      end if
     end if
      
    case(2)
     if(g_cell%edge(num).eq.1) then
      direction='right  '
     else
      if(g_cell%edge(num).eq.2) then
       direction='down  '
      else
       call error_message
      end if
     end if   
     
    case(3)
     if(g_cell%edge(num).eq.2) then
      direction='up     '
     else
      if(g_cell%edge(num).eq.3) then
       direction='right '
      else
       call error_message
      end if
     end if   
    
    case(4)
     if(g_cell%edge(num).eq.3) then
      direction='left   '
     else
      if(g_cell%edge(num).eq.4) then
       direction='up    '
      else
       call error_message
      end if
     end if   
    
   end select
    
  else
   
   select case(single)
    
    case(1)
     if(g_cell%edge(num).eq.1) then
      direction='right '
     else
      if(g_cell%edge(num).eq.4) then
       direction='up    '
      else
       call error_message
      end if
     end if
      
    case(2)
     if(g_cell%edge(num).eq.1) then
      direction='left   '
     else
      if(g_cell%edge(num).eq.2) then
       direction='up    '
      else
       call error_message
      end if
     end if   
     
    case(3)
     if(g_cell%edge(num).eq.2) then
      direction='down  '
     else
      if(g_cell%edge(num).eq.3) then
       direction='left  '
      else
       call error_message
      end if
     end if   
    
    case(4)
     if(g_cell%edge(num).eq.3) then
      direction='right '
     else
      if(g_cell%edge(num).eq.4) then
       direction='down  '
      else
       call error_message
      end if
     end if
    	 
   end select
    
  end if

end select

! Finally, determine whether the two positions crossing.
select case(direction)
 case('left  ','down  ')
  if(g_celln%dis_posi(numn).gt.g_cell%dis_posi(num)) cross=.true.
 case('right ','up    ')
  if(g_celln%dis_posi(numn).lt.g_cell%dis_posi(num)) cross=.true.   
 case default; call error_message
end select 
    
end subroutine if_crossing


end subroutine adjust_crossing


subroutine modify_neigh(tempx,numx)

implicit none

type(auxi_crit_cell), pointer :: tempx
integer, intent(in) ::numx

type(auxi_crit_cell), pointer :: tempxx
integer :: i, j, i1, j1
real(8) :: difference
  
i=tempx%g_cell%x_idx; j=tempx%g_cell%y_idx
  
if(tempx%partner.ne.'single  ') then 
   
 if(tempx%partner.eq.'previous') then
  tempxx=>tempx%previous
 else
  tempxx=>tempx%next
 end if
  
 i1=tempxx%g_cell%x_idx; j1=tempxx%g_cell%y_idx 
 if(tempxx%g_cell%g_type.eq.'xx ') then
  difference=dfloat(i-i1)
 else
  if(tempxx%g_cell%g_type.eq.'yy ') then
   difference=dfloat(j-j1)
  else
   call error_message
  end if
 end if
  
 tempxx%g_cell%dis_posi(numx)=tempx%g_cell%dis_posi(numx)+difference
  
end if
   
select case(numx)
   
 case(1)
  
  if(tempx%partner.eq.'previous') tempx=>tempx%previous
  
  tempx=>tempx%previous
  tempx%g_cell%dis_posi(2)=tempx%next%g_cell%dis_posi(1)
   
  i=tempx%g_cell%x_idx; j=tempx%g_cell%y_idx
  
  if(tempx%partner.eq.'next    ') call error_message
   
  if(tempx%partner.eq.'previous') then   
   tempxx=>tempx%previous  
   i1=tempxx%g_cell%x_idx; j1=tempxx%g_cell%y_idx 
   if(tempxx%g_cell%g_type.eq.'xx ') then
    difference=dfloat(i-i1)
   else
    if(tempxx%g_cell%g_type.eq.'yy ') then
     difference=dfloat(j-j1)
    else
     call error_message
    end if
   end if
   tempxx%g_cell%dis_posi(2)=tempx%g_cell%dis_posi(2)+difference
  end if
   
 case(2)
  
  if(tempx%partner.eq.'next    ') tempx=>tempx%next
  
  tempx=>tempx%next
  tempx%g_cell%dis_posi(1)=tempx%previous%g_cell%dis_posi(2)
  
  i=tempx%g_cell%x_idx; j=tempx%g_cell%y_idx
  
  if(tempx%partner.eq.'previous') call error_message
   
  if(tempx%partner.eq.'next    ') then   
   tempxx=>tempx%next  
   i1=tempxx%g_cell%x_idx; j1=tempxx%g_cell%y_idx 
   if(tempxx%g_cell%g_type.eq.'xx ') then
    difference=dfloat(i-i1)
   else
    if(tempxx%g_cell%g_type.eq.'yy ') then
     difference=dfloat(j-j1)
    else
     call error_message
    end if
   end if
   tempxx%g_cell%dis_posi(1)=tempx%g_cell%dis_posi(1)+difference
  end if
   
 case default; call error_message
  
end select

end subroutine modify_neigh
   

end module manipulation_2



