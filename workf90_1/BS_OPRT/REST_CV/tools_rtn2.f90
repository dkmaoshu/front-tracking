module tools_for_reset
! This module provides the tools needed in the reset of 
! discontinuity curves.

use solu_comput
! 'solu_comp.f90'

use tools
! 'tools.f90'

implicit none
character*6, dimension(2) :: move_dir
! Tells the move direction of the two discontinuity positions, 
! either from left to right 'l_to_r' or from right to left
! 'r_to_l'.

character*6 :: happen
! Tells whether the critical cell is going to move, split, or
! remain its status.

public  move_dir, happen, set_move_1, geo_split, geo_move
private adjust_posi

		  
contains


subroutine set_move_1(temp_set,moved,happen)
! Determine whether the two discontinuity positions in the 
! auxiliary critical cell should move and in which direction 
! according to their positions.

implicit none
type(auxi_crit_cell), pointer :: temp_set
! The pointer pointing to the auxiliary critical cell under 
! concern.
logical, dimension(2), intent(out)  :: moved
! Tell whether the two discontinuity positions will move.
character*6, intent(out) :: happen

type(auxi_crit_cell), pointer :: temp_stack, temp_partner
type(geo_info) :: g_cell, g_stack, g_partner
! The geometrical cell of the critical cell under concern the 
! left or right stacked critical cell, and the partner critical
! cell.
logical, dimension(2) :: be_out
! Tell whether the position should be out of the critical cell. This information
! is used for partnered critical cells
character*6, dimension(2) :: out_dir
! Tell the oriental should-be-out directions.
integer :: i, j, sgl_c, i0, j0, i1, j1
character*8 :: partner, direction
logical :: crossed, if_modify
integer :: i_present, j_present, i_stacked, j_stacked, edge_between


move_dir=' '; g_cell=temp_set%g_cell; moved=.false.
partner=temp_set%partner; nullify(temp_partner)
call clean_up_g_cell(g_partner)
be_out=.false.; out_dir='      '

! Determine if-should-be-out information for the two discontinuity positions.
i0=g_cell%x_idx; j0=g_cell%y_idx;
i1=error_index; j1=error_index
select case(partner)
 case('previous')
  temp_partner=>temp_set%previous
  g_partner=temp_partner%g_cell; be_out(1)=.true.
  i1=g_partner%x_idx; j1=g_partner%y_idx
 case('next    ')
  temp_partner=>temp_set%next
  g_partner=temp_partner%g_cell; be_out(2)=.true.
  i1=g_partner%x_idx; j1=g_partner%y_idx
end select 

do i=1,2
 if(be_out(i)) then
  if(i1+j1.gt.i0+j0) then
   if(g_cell%point(1).eq.'left  ') then
    out_dir(i)='l_to_r'
   else
    out_dir(i)='r_to_l'
   end if
  else
   if(g_cell%point(1).eq.'left  ') then
    out_dir(i)='r_to_l'
   else
    out_dir(i)='l_to_r'
   end if
  end if
 end if
end do

! First, determine whether the two discontinity positions will 
! move or stay according to their locations.
do i=1,2
 if(dabs(g_cell%dis_posi(i)).gt.0.5d0) then
  select case(g_cell%edge(i))
   case(1,2)
    if(g_cell%dis_posi(i).gt.0.5d0) then
     if(g_cell%point(g_cell%edge(i)).eq.'left  ') then
      move_dir(i)='l_to_r'
     else
      move_dir(i)='r_to_l'
     end if
    else
	 if(g_cell%point(g_cell%edge(i)).eq.'left  ') then
      move_dir(i)='r_to_l'
     else
      move_dir(i)='l_to_r'
     end if
    end if
   case(3,4)
    if(g_cell%dis_posi(i).gt.0.5d0) then
     if(g_cell%point(g_cell%edge(i)).eq.'left  ') then
      move_dir(i)='r_to_l'
     else
      move_dir(i)='l_to_r'
     end if
    else
     if(g_cell%point(g_cell%edge(i)).eq.'left  ') then
      move_dir(i)='l_to_r'
     else
      move_dir(i)='r_to_l'
     end if
    end if
  end select
 end if
end do

do i=1,2
 if(.not.be_out(i)) then
  if(move_dir(i).ne.'      ') moved(i)=.true.
 else
  if(move_dir(i).eq.'l_to_r'.and.out_dir(i).eq.'r_to_l'.or. &
     move_dir(i).eq.'r_to_l'.and.out_dir(i).eq.'l_to_r')  &
   call error_message
  if(move_dir(i).eq.'      ') moved(i)=.true.
 end if
end do

! Modify the move-split-stay situation based on the stacked
! situation.
do i=1,2
 nullify(temp_stack); call clean_up_g_cell(g_stack)
 call determine_index_from_partner(temp_set,i,i_present,j_present)
 call index_n_move_direction(i,i_present,j_present,direction)
 if(.not.be_out(i)) then
  select case(move_dir(i))
   case('l_to_r')
    if(temp_set%r_stk%cv_nmb.gt.0) &
     call visit(temp_set%r_stk,temp_stack)
   case('r_to_l')
    if(temp_set%l_stk%cv_nmb.gt.0) &
     call visit(temp_set%l_stk,temp_stack)
  end select
 else
  if(be_out(i).and.move_dir(i).eq.'      ') then
   if(out_dir(i).eq.'l_to_r') then
    if(temp_partner%l_stk%cv_nmb.gt.0)  &
     call visit(temp_partner%l_stk,temp_stack)
   else
    if(temp_partner%r_stk%cv_nmb.gt.0)  &
     call visit(temp_partner%r_stk,temp_stack)
   end if
  end if
 end if
  
 if(associated(temp_stack)) then
  g_stack=temp_stack%g_cell
   
  if(g_stack%g_type.eq.g_cell%g_type.or.temp_stack%partner.eq.'single  ') then    
   
   do j=1,2
    call determine_index_from_partner(temp_stack,j,i_stacked,j_stacked)
    if_modify=.false.
    
    if(g_stack%edge(j).eq.g_cell%edge(i)) then
     select case(g_cell%edge(i))
      case(1,3)
       if(direction.eq.'positive'.and.i_present.gt.i_stacked.or. &
          direction.eq.'negative'.and.i_present.lt.i_stacked) then
        if_modify=.true.; exit
       end if       
      case(2,4)
       if(direction.eq.'positive'.and.j_present.gt.j_stacked.or. &
          direction.eq.'negative'.and.j_present.lt.j_stacked) then
        if_modify=.true.; exit
       end if
     end select 
    end if
   end do
  else
   call find_edge_in_between_4_partnered(temp_stack,edge_between)
   if(edge_between.eq.g_cell%edge(i)) then
    select case(g_cell%edge(i))
     case(1,3)
      if(direction.eq.'positive'.and.i_present.gt.g_stack%x_idx.or. &
         direction.eq.'negative'.and.i_present.lt.g_stack%x_idx) if_modify=.true.
     case(2,4)
      if(direction.eq.'positive'.and.j_present.gt.g_stack%y_idx.or. &
         direction.eq.'negative'.and.j_present.lt.g_stack%y_idx) if_modify=.true.
    end select     
   end if
  end if
    
  if(if_modify) then 
   if(.not.be_out(i)) then
    move_dir(i)='      '; moved(i)=.false.
    g_cell%dis_posi(i)=dsign(0.4999999d0,g_cell%dis_posi(i))
   else
    move_dir(i)=out_dir(i); moved(i)=.false.
    g_cell%dis_posi(i)=dsign(0.5000001d0,g_cell%dis_posi(i))
   end if
  end if
    
 end if
end do

! Fix unconsistency in movement of 'xy'-type critical cells.
if(g_cell%g_type.eq.'xy '.and.move_dir(1).ne.move_dir(2)) then
 call find_single(g_cell,sgl_c)
 crossed=.true.
 do i=1,2
  if(move_dir(i).eq.'l_to_r'.and.g_cell%point(sgl_c).eq.'right '.or. &
     move_dir(i).eq.'r_to_l'.and.g_cell%point(sgl_c).eq.'left  ') then
   crossed=.false.; exit
  end if
 end do
 if(.not.crossed) then
  move_dir='      '; moved=.false.
  call adjust_posi(g_cell,sgl_c,1,'no ')
  call adjust_posi(g_cell,sgl_c,2,'no ')
  temp_set%g_cell=g_cell
 end if
end if

! The following is to determine what is going to happen to the 
! critical cell under concern.
select case(g_cell%g_type)
 case('xx ', 'yy ')
  if(move_dir(1).eq.move_dir(2)) then
   if(move_dir(1).eq.'      ') then
    happen='remain'
   else
    happen='move  '
   end if
  else
   happen='split '
  end if
 case('xy ')
  if(move_dir(1).ne.move_dir(2)) then
   happen='split '
  else
   call find_single(g_cell,sgl_c)
   if(move_dir(1).eq.'l_to_r'.and.g_cell%point(sgl_c).eq.'right '.or. &
      move_dir(1).eq.'r_to_l'.and.g_cell%point(sgl_c).eq.'left  ') then
    happen='move  '
   else
    if(move_dir(1).eq.'      ') then
     happen='remain'
    else
     happen='split '
    end if
   end if
  end if		  
end select

temp_set%g_cell=g_cell


contains


subroutine determine_index_from_partner(temp,edge_num,i,j)

implicit none
type(auxi_crit_cell) :: temp
! The present auxiliary critical cell.
integer, intent(in) :: edge_num
integer, intent(out) :: i, j

type(geo_info) :: g_present, g_partner
logical :: if_out

i=error_index; j=error_index
g_present=temp%g_cell; call clean_up_g_cell(g_partner)
if_out=.false.

if(temp%partner.eq.'previous'.and.edge_num.eq.1) if_out=.true.
if(temp%partner.eq.'next    '.and.edge_num.eq.2) if_out=.true.
if(.not.if_out) then
 select case(g_present%edge(edge_num))
  case(1,3); i=g_present%x_idx
  case(2,4); j=g_present%y_idx
 end select
else
 select case(temp%partner)
  case('previous'); g_partner=temp%previous%g_cell
  case('next    '); g_partner=temp%next%g_cell
  case default; call error_message; pause
 end select
 select case(g_present%edge(edge_num))
  case(1,3); i=g_partner%x_idx
  case(2,4); j=g_partner%y_idx
 end select
end if

end subroutine determine_index_from_partner


subroutine index_n_move_direction(edge_num,i,j,direction)

implicit none
integer, intent(in) :: edge_num
integer, intent(inout) :: i, j
character*8, intent(out) :: direction

type(geo_info) :: g_present
integer :: ii, jj

ii=error_index; jj=error_index; direction='        '
g_present=temp_set%g_cell

select case(g_present%edge(edge_num))
 case(1)
  if(move_dir(edge_num).eq.'l_to_r'.and.g_present%point(1).eq.'left  ' &
     .or.move_dir(edge_num).eq.'r_to_l'.and.g_present%point(1).eq. &
	 'right ') then
   ii=g_present%x_idx+1
  else
   ii=g_present%x_idx-1
  end if
  if(move_dir(edge_num).eq.'      ') ii=g_present%x_idx
 case(2)
  if(move_dir(edge_num).eq.'l_to_r'.and.g_present%point(2).eq.'left  ' &
     .or.move_dir(edge_num).eq.'r_to_l'.and.g_present%point(2).eq. &
	 'right ') then
   jj=g_present%y_idx+1
  else
   jj=g_present%y_idx-1
  end if
  if(move_dir(edge_num).eq.'      ') jj=g_present%y_idx
 case(3)
  if(move_dir(edge_num).eq.'l_to_r'.and.g_present%point(4).eq.'left  ' &
     .or.move_dir(edge_num).eq.'r_to_l'.and.g_present%point(4).eq. &
	 'right ') then
   ii=g_present%x_idx+1
  else
   ii=g_present%x_idx-1
  end if
  if(move_dir(edge_num).eq.'      ') ii=g_present%x_idx
 case(4)
  if(move_dir(edge_num).eq.'l_to_r'.and.g_present%point(1).eq.'left  ' &
     .or.move_dir(edge_num).eq.'r_to_l'.and.g_present%point(1).eq. &
	 'right ') then
   jj=g_present%y_idx+1
  else
   jj=g_present%y_idx-1
  end if
  if(move_dir(edge_num).eq.'      ') jj=g_present%x_idx
end select 

select case(g_present%edge(edge_num))
 case(1,3)
  if(ii.gt.i) direction='positive'
  if(ii.lt.i) direction='negative'
 case(2,4)
  if(jj.gt.j) direction='positive'
  if(jj.lt.j) direction='negative'
end select 
  
i=ii; j=jj

end subroutine index_n_move_direction


end subroutine set_move_1 


subroutine adjust_posi(g_cell,sgl_c,number,crossed)

implicit none
type(geo_info), intent(inout) :: g_cell
integer, intent(in) :: sgl_c, number
character*3, intent(in) :: crossed

select case(crossed)
 
 case('yes')

  select case(sgl_c)
   case(1)
  	if(g_cell%dis_posi(number).ge.-0.5d0) then
     g_cell%dis_posi(number)=-0.5000001d0
    end if
   case(2)
    if(g_cell%edge(number).eq.1.and.g_cell%dis_posi(number).le.0.5d0) then
     g_cell%dis_posi(number)=0.5000001d0
    end if
    if(g_cell%edge(number).eq.2.and.g_cell%dis_posi(number).ge.-0.5d0) then
     g_cell%dis_posi(number)=-0.5000001d0
    end if
   case(3)
  	if(g_cell%dis_posi(number).le.0.5d0) then
     g_cell%dis_posi(number)=0.5000001d0
    end if
   case(4)
    if(g_cell%edge(number).eq.4.and.g_cell%dis_posi(number).le.0.5d0) then
     g_cell%dis_posi(number)=0.5000001d0
    end if
    if(g_cell%edge(number).eq.3.and.g_cell%dis_posi(number).ge.-0.5d0) then
     g_cell%dis_posi(number)=-0.5000001d0
    end if
  end select

 case('no ')

  select case(sgl_c)
   case(1)
  	if(g_cell%dis_posi(number).le.-0.5d0) then
     g_cell%dis_posi(number)=-0.4999999d0
    end if
   case(2)
    if(g_cell%edge(number).eq.1.and.g_cell%dis_posi(number).ge.0.5d0) then
     g_cell%dis_posi(number)=0.4999999d0
    end if
    if(g_cell%edge(number).eq.2.and.g_cell%dis_posi(number).le.-0.5d0) then
     g_cell%dis_posi(number)=-0.4999999d0
    end if
   case(3)
  	if(g_cell%dis_posi(number).ge.0.5d0) then
     g_cell%dis_posi(number)=0.4999999d0
    end if
   case(4)
    if(g_cell%edge(number).eq.4.and.g_cell%dis_posi(number).ge.0.5d0) then
     g_cell%dis_posi(number)=0.4999999d0
    end if
    if(g_cell%edge(number).eq.3.and.g_cell%dis_posi(number).le.-0.5d0) then
     g_cell%dis_posi(number)=-0.4999999d0
    end if
  end select

end select

end subroutine adjust_posi


subroutine geo_split(stencil,g_cell,g_new,num,curve_coefficients)
! Geometrically split the auxiliary critical cell 'g_cell' 
! into two auxiliary critical cells 'g_cell' and 'g_new' when 
! the 'num'th discontinuity position goes out of the cell.

implicit none

real(8), dimension(3,2), intent(in) :: stencil
type(geo_info), intent(inout) :: g_cell
type(geo_info), intent(out) :: g_new
integer, intent(in) :: num
real(8), dimension(3), intent(out) :: curve_coefficients

real(8), dimension(2) :: pt
integer :: edge, replace, edge_o, i
integer :: pt_num, edge_num, size
real(8) :: value1, value2, dis_posi_o, dif_idx, root, right_term
real(8), dimension(:), allocatable :: values, points, coefficients
real(8), dimension(:), allocatable :: roots

call clean_up_g_cell(g_new)
edge_num=g_cell%edge(num)

! Choose the points corresponding to the discontinuity curve 
! segment and the cell-edge the segment crosses.
if(stencil(3,1).gt.0.9d0*error_data) then
 size=3
else
 size=2
end if
allocate(values(size),points(size),coefficients(size),roots(size-1))
values=error_data; points=error_data; coefficients=error_data; roots=error_data
select case(edge_num)
 case(1,3)
  values=stencil(1:size,1); points=stencil(1:size,2)
 case(2,4)
  values=stencil(1:size,2); points=stencil(1:size,1)
 case default; call error_message
end select
select case(move_dir(num))
 case('l_to_r')
  if(g_cell%point(edge_num).eq.'left  ') then
   replace=1
  else
   replace=-1
  end if
 case('r_to_l')
  if(g_cell%point(edge_num).eq.'left  ') then
   replace=-1
  else
   replace=1
  end if
end select 
edge=cycle_4(edge_num+replace)

! Find the intersection point.
call polynomial(values,points,coefficients)
select case(edge)
 case(1,4); right_term=-0.5d0
 case(2,3); right_term=0.5d0
end select
call polynomial_roots(coefficients,right_term,roots)
select case(size-1)
 case(1)
  select case(edge)
   case(2,4); pt(1)=right_term; pt(2)=roots(1)
   case(1,3); pt(1)=roots(1); pt(2)=right_term
  end select
 case(2)
  if(dabs(roots(1)).gt.dabs(roots(2))) then
   root=roots(2)
  else
   root=roots(1)
  end if
  select case(edge)
   case(2,4); pt(1)=right_term; pt(2)=root
   case(1,3); pt(1)=root; pt(2)=right_term
  end select
end select
! Restrict the intersection point when one of the components of
! 'move_dir' describes a fake move.

select case(edge_num)
 case(1)
  pt_num=2; value1=-0.4999999d0; value2=0.4999999d0
 case(2)
  pt_num=1; value1=0.4999999d0; value2=-0.4999999d0
 case(3)
  pt_num=2; value1=0.4999999d0; value2=-0.4999999d0
 case(4)
  pt_num=1; value1=-0.4999999d0; value2=0.4999999d0
end select
if(dabs(g_cell%dis_posi(num)).le.0.5d0) then
 if(pt(pt_num).gt.0.5d0) pt(pt_num)=value1
 if(pt(pt_num).gt.0.5d0) pt(pt_num)=value2
end if

! Update the original cell 'g_cell'.
edge_o=g_cell%edge(num); g_cell%edge(num)=edge
dis_posi_o=g_cell%dis_posi(num)
select case(edge)
 case(1,3); g_cell%dis_posi(num)=pt(1)
 case(2,4); g_cell%dis_posi(num)=pt(2)
end select
call find_type_from_edges(g_cell)
call find_points_from_edges(g_cell)

! Define the new cell 'g_new'.
g_new%edge(num)=edge_o
g_new%edge(3-num)=neighbor_edge(edge)
g_new%x_idx=g_cell%x_idx; g_new%y_idx=g_cell%y_idx
select case(edge)
 case(1); g_new%y_idx=g_new%y_idx-1
 case(2); g_new%x_idx=g_new%x_idx+1
 case(3); g_new%y_idx=g_new%y_idx+1
 case(4); g_new%x_idx=g_new%x_idx-1
end select
dif_idx=g_new%x_idx+g_new%y_idx-g_cell%x_idx-g_cell%y_idx
g_new%dis_posi(num)=dis_posi_o-dif_idx
g_new%dis_posi(3-num)=g_cell%dis_posi(num)
call find_type_from_edges(g_new)
call find_points_from_edges(g_new)

curve_coefficients=error_data
do i=1, size
 curve_coefficients(i)=coefficients(i)
end do

deallocate(values,points,coefficients,roots)

end subroutine geo_split


subroutine geo_move(g_cell)
! Geometrically move the auxiliary critical cell when its two
! discontinuity positions go out of the cell in the same direction.

implicit none
type(geo_info), intent(inout) :: g_cell

type(geo_info) :: gg
integer :: single, i

gg=g_cell
select case(g_cell%g_type)

 case('xx','yy')
  if(g_cell%dis_posi(1).gt.0.5d0) then
   gg%dis_posi=gg%dis_posi-1.0d0
   if(g_cell%g_type.eq.'xx') then
    gg%x_idx=gg%x_idx+1
   else
    gg%y_idx=gg%y_idx+1
   end if
  else
   gg%dis_posi=gg%dis_posi+1.0d0
   if(g_cell%g_type.eq.'xx') then
    gg%x_idx=gg%x_idx-1
   else
    gg%y_idx=gg%y_idx-1
   end if
  end if
  
 case('xy')
  call find_single(g_cell,single)
  select case(single)
   case(1) 
    gg%dis_posi=g_cell%dis_posi+1.0d0
    gg%x_idx=gg%x_idx-1; gg%y_idx=gg%y_idx-1
   case(2)
    do i=1,2    	 	  	  
     if(g_cell%edge(i).eq.1) gg%dis_posi(i)=g_cell%dis_posi(i)-1.d0
     if(g_cell%edge(i).eq.2) gg%dis_posi(i)=g_cell%dis_posi(i)+1.d0
    end do
    gg%x_idx=gg%x_idx+1; gg%y_idx=gg%y_idx-1
   case(3)
    gg%dis_posi=g_cell%dis_posi-1.0d0
    gg%x_idx=gg%x_idx+1; gg%y_idx=gg%y_idx+1
   case(4)
    do i=1,2    	 	  	  
     if(g_cell%edge(i).eq.3) gg%dis_posi(i)=g_cell%dis_posi(i)+1.d0
     if(g_cell%edge(i).eq.4) gg%dis_posi(i)=g_cell%dis_posi(i)-1.d0
    end do
    gg%x_idx=gg%x_idx-1; gg%y_idx=gg%y_idx+1
  end select
  do i=1,2
   gg%edge(i)=neighbor_edge(g_cell%edge(i))
  end do
   
  call find_points_from_edges(gg)

!%%%%%%%%%%%%%%%%% patch %%%%%%%%%%%%
! The points' side of a moved $xy$-type discontinyity cell should be
! reversed.
  do i=1,4
   if(gg%point(i).eq.'left  ') then
    gg%point(i)='right '
   else
    if(gg%point(i).eq.'right ') then
	 gg%point(i)='left  '
    else
     call error_message
    end if
   end if
  end do    
!%%%%%%%%%%%%%%%%% patch %%%%%%%%%%%%
  
end select

g_cell=gg

end subroutine geo_move

end module tools_for_reset