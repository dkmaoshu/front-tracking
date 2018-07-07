module pick_discontinuity_positions

use auxiliary_discontinuity_curves
! 'auxi_cv.f90'

implicit none

type(auxi_crit_cell), pointer :: tempa
type(geo_info) :: g_cell, g_current

public  pick_positions, tempa
private g_cell


contains


subroutine pick_positions(starting_posi,positions,point_number, &
                          direction)

implicit none
integer, intent(in) :: starting_posi
! The number of the starting position.
real(8), dimension(3,2), intent(out) :: positions
! The picked the discontinuity positions.
integer, intent(out) :: point_number
! The number of the points being picked.
character*8, intent(out) :: direction
! The direction the picking proceeds.

type(auxi_crit_cell), pointer :: tempa_current
real(8), dimension(2) :: pt
integer :: picking_number
! The number of the discontinuity position going to be picked.
logical :: if_end, if_crossed, if_length_ok, if_other_end

positions=error_data
g_cell=tempa%g_cell; direction='        '
call clean_up_g_cell(g_current); if_end=.false.

! Pick the first position.
call pick_point(g_cell,starting_posi,pt)
point_number=1; positions(1,:)=pt

! The following determines the picking direction and so on.
select case(starting_posi)
 
 case(1)
  select case(g_cell%wind(1))
   case('in ')
    if(associated(tempa%p_nxt_dr)) then
     tempa_current=>tempa%p_nxt_dr
     direction='previous'
     picking_number=1
    else
     tempa_current=>tempa
     direction='next    '
     picking_number=2
     if_end=.true.
    end if
      
   case('out')
    tempa_current=>tempa
    direction='next    '
    picking_number=2

   case default; call error_message
  end select
 
 case(2)
  select case(g_cell%wind(2))
   case('in ')
    if(associated(tempa%n_nxt_dr)) then
     tempa_current=>tempa%n_nxt_dr
     direction='next    '
     picking_number=2
    else
     tempa_current=>tempa
     direction='previous'
     picking_number=1
     if_end=.true.
    end if
      
   case('out')
    tempa_current=>tempa
    direction='previous'
    picking_number=1

   case default; call error_message
  end select

 case default; call error_message

end select

do while(point_number.lt.3)
 
 g_current=tempa_current%g_cell

! The following deals with critical cells of 'xy '-type either 
! whose discontinuity positions have crossed the corner or whose
! discontinuity segment is too short.
 if(g_current%g_type.eq.'xy ') then
  call check_crossed(g_current,if_crossed)
  call check_length(g_current,if_length_ok)

  if(if_crossed.or..not.if_length_ok) then
   select case(picking_number)

    case(1)
     if(associated(tempa_current%p_nxt_dr))	then
      tempa_current=>tempa_current%p_nxt_dr
     else
      if(if_end) then
       return
      else
       select case(starting_posi)
        case(1); tempa_current=>tempa
        case(2)
         if(associated(tempa%n_nxt_dr)) then
          tempa_current=>tempa%n_nxt_dr
    	 else
          return
         end if
       end select
       picking_number=2
       if(point_number.eq.1) then
        direction='next    '
       else
        direction='bet_prev'
       end if
      end if
     end if

    case(2)
     if(associated(tempa_current%n_nxt_dr))	then
      tempa_current=>tempa_current%n_nxt_dr
     else
      if(if_end) then
       return
      else
       select case(starting_posi)
        case(2); tempa_current=>tempa
        case(1)
         if(associated(tempa%p_nxt_dr)) then
          tempa_current=>tempa%p_nxt_dr
         else
          return
         end if
       end select
       picking_number=1
       if(point_number.eq.1) then
        direction='previous'
       else
        direction='bet_next'
       end if
      end if
     end if

    case default; call error_message

   end select

   g_current=tempa_current%g_cell
  end if

 end if

! Picks a position.
 call pick_point(g_current,picking_number,pt)
 point_number=point_number+1
 pt(1)=pt(1)+dfloat(g_current%x_idx-g_cell%x_idx)
 pt(2)=pt(2)+dfloat(g_current%y_idx-g_cell%y_idx)
 positions(point_number,:)=pt

! Go to the next critical cell.
 if(point_number.eq.3) return
  
 if_other_end=.false.
 select case(picking_number)
 
  case(1)
   if(associated(tempa_current%p_nxt_dr)) then
    tempa_current=>tempa_current%p_nxt_dr
   else
    if_other_end=.true.
   end if

  case(2)
   if(associated(tempa_current%n_nxt_dr)) then
    tempa_current=>tempa_current%n_nxt_dr
   else
    if_other_end=.true.
   end if

  case default; call error_message
 end select
  
 if(if_other_end) then   
  if(if_end) then
   return
  else
   select case(picking_number)

    case(1)
     select case(starting_posi)
      case(1); tempa_current=>tempa
      case(2)
       if(associated(tempa%n_nxt_dr)) then
        tempa_current=>tempa%n_nxt_dr
       else
        return
       end if
     end select
     direction='bet_prev'
     picking_number=2
        	  
    case(2)	  	  	    
     select case(starting_posi)
      case(2); tempa_current=>tempa
      case(1)
       if(associated(tempa%p_nxt_dr)) then
        tempa_current=>tempa%p_nxt_dr
       else
        return
       end if
     end select
     direction='bet_next'
     picking_number=1	  	  				 				 			   	  
       
   end select
  end if
 end if

end do


contains


subroutine check_crossed(g_cell,if_crossed)

implicit none
type(geo_info), intent(in) :: g_cell
logical, intent(out) :: if_crossed

integer :: single, i

if(g_cell%g_type.ne.'xy ') call error_message
if_crossed=.false.

call find_single(g_cell,single)

select case(single)
 
 case(1)
  do i=1,2
   if(g_cell%dis_posi(i).lt.-0.5d0) then
    if_crossed=.true.; return
   end if
  end do
 
 case(2)
  do i=1,2
   if(g_cell%edge(i).eq.1.and.g_cell%dis_posi(i).gt.0.5d0.or. &
      g_cell%edge(i).eq.2.and.g_cell%dis_posi(i).lt.-0.5d0) then
    if_crossed=.true.; return
   end if
  end do

 case(3)
  do i=1,2
   if(g_cell%dis_posi(i).gt.0.5d0) then
    if_crossed=.true.; return
   end if
  end do
 
 case(4)
  do i=1,2
   if(g_cell%edge(i).eq.4.and.g_cell%dis_posi(i).gt.0.5d0.or. &
      g_cell%edge(i).eq.3.and.g_cell%dis_posi(i).lt.-0.5d0) then
    if_crossed=.true.; return
   end if
  end do

end select

end subroutine check_crossed


subroutine check_length(g_cell,if_ok)

implicit none
type(geo_info), intent(in) :: g_cell
logical, intent(out) :: if_ok

real(8), dimension(2) :: pt1, pt2
real(8) :: length

if_ok=.true.

call pick_point(g_cell,1,pt1)
call pick_point(g_cell,2,pt2)
length=length_of_segment(pt1,pt2)

if(length.lt.0.2d0) if_ok=.false.

end subroutine check_length


end subroutine pick_positions


end module pick_discontinuity_positions