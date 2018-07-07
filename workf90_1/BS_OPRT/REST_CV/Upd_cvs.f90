module update_curves
! This module is for update of the auxiliary discontinuity curve
! by moving and spliting critical cells if necessary.

use upd_cv_cases_1
! 'updcvcas_1.f90'

use remove_tangles
! 'rmv_tgls.f90'

implicit none

integer :: head_mark, end_mark
logical :: head_switch

public  update_acvs
private update_acv, set_move_2, move_update, move_dir, happen, &
        move, split, side_neighbor, head_mark, head_switch, &
		ulx_exd_p, uly_exd_p, urx_exd_p, ury_exd_p, fix_heads, &
		fix_remaining_tangled


contains


subroutine update_acvs
! Update all auxiliary discontinuity curves.

implicit none
integer :: i, start_index, end_index, step

logical :: possible_remaining

! type(grid_cell) :: rg

! call check_list_ac(1,'down  ')
! rg=ggd_cell(15,-6)

select case(order_of_curve_reset)
 case('right_way   '); start_index=1; end_index=cvn; step=1
 case('inverse_way '); start_index=cvn; end_index=1; step=-1
 case default; call error_message
end select

do i=start_index, end_index, step
!do i=1,cvn
 if(acvv(i)%status.eq.'awake ') then

  if(acvv(i)%cv_type.eq.'circular') call head_shift(acvv(i))
! The make sure that the two heads of the curve are at a right place
! so that there will not be complications caused in the reset.

!   call check_list_ac(1,'down  ')

  call update_acv(acvv(i))

!   call check_list_ac(1,'down  ')
!   rg=ggd_cell(11,-6)
 
  if(acvv(i)%cv_type.eq.'circular') call fix_heads(acvv(i)) 

!   call check_list_ac(1,'down  ')

  possible_remaining=.true.
  do while(possible_remaining)
   call fix_remaining_tangled(acvv(i),possible_remaining)
  end do

!   call check_list_ac(1,'down  ')

 end if
end do

end subroutine update_acvs


subroutine update_acv(acv)
! Update a single auxiliary discontinuity curve.

implicit none
type(auxi_discv), intent(inout) :: acv
! The single auxiliary discontinuity curve to be updated.

integer :: i, ii
type(geo_info) :: gn_cell
type(state), dimension(-1:1) :: ulepx, ulepy, urepx, urepy

type(auxi_crit_cell), pointer :: boundary_end
! integer ::idxx

! integer, dimension(3) :: stencil_disp
   
! character*6 :: rg

! type(auxi_crit_cell), pointer :: tt
! type(adss_info) :: adss, adss_l, adss_r

! rg=ggd_cell(-2,-7)%region

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

do ii=-1, 1
 ulx_exd_p(ii)=error_data; uly_exd_p(ii)=error_data
 urx_exd_p(ii)=error_data; ury_exd_p(ii)=error_data
end do

do while(associated(temp).and.head_switch) 
 g_cell=temp%g_cell; p_cell=temp%p_cell
 ulepx=temp%ulx_exd; ulepy=temp%uly_exd
 urepx=temp%urx_exd; urepy=temp%ury_exd
  
 call clean_up_g_cell(gn_cell)
 if(associated(temp%next)) gn_cell=temp%next%g_cell
  
!  rg=ggd_cell(-2,-7)%region
! Determine the move situation.
  
!  call check_list_ac(1,'down  ')
  
!  if(associated(temp%n_nxt_dr)) then
!   tt=>temp%n_nxt_dr
!   ggg_cell=tt%g_cell
!   stencil_disp=tt%stencil_disp
!  end if
  
 call set_move_2(temp)
  
!  if(associated(temp%previous)) ggg_cell=temp%previous%g_cell

!  call check_list_ac(1,'up    ')
   
!  Update the auxiliary critical cell under concern.
 select case(happen)
  case('move  '); call move_update

!   call check_list_ac(1,'up    ')  
!   if(g_cell%x_idx.eq.3) print*, ggd_cell(3,0)%region
  
  case('split ')
   if(temp%partner.ne.'next    ') then 
    do i=1,2
     if(move_dir(i).ne.'      ') call split_update(i)
    end do    
   else
    do i=2,1,-1
     if(move_dir(i).ne.'      ') call split_update(i)
    end do    
   end if
   
!    call check_list_ac(1,'up    ')
  
  case('remain') 
   temp%updated='yes'; temp%partner='single  '
 end select

 temp%p_cell%or_state%gamma=error_data
  
!  call check_list_ac(1,'down  ')
!  call visit(adss_info(1,1),tt)
!  adss=temp%address; adss_l=temp%l_stk; adss_r=temp%r_stk
  
 ulx_exd_p=ulepx; uly_exd_p=ulepy
 urx_exd_p=urepx; ury_exd_p=urepy
   
 temp=>temp%next
 if(associated(temp)) then
  
!   idxx=temp%address%idx
  
  head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
  
end do

if(acv%cv_type.eq.'circular') then
 acv%begin_boundary%next=>acv%begin%next
 acv%eend_boundary%previous=>acv%eend%previous
end if

end subroutine update_acv


subroutine move_update
! This subroutine updates the auxiliary critical cell that should
! move and its neighbor-relation and grid-map around the cell.

implicit none

type(state) :: u_new
! The data of the smooth solution left behind in the old cell
! after the move.
integer, dimension(2) :: cross_o, cross_n
! The corner points the move will cross, 'cross_o', viewed in
! the old cell and 'cross_n', viewed in the new cell. If the
! moved auxiliary critical cell is of either 'xx' or 'yy' type,
! it corsses two corner points, otherwise it crosses only one
! corner point.
type(auxi_crit_cell), pointer :: temp_g, temp_partner, temp_gg
! Pointer to the spatial neighboring auxiliary critical cell 
! in the grid cell.
character*6 :: this_side, that_side, head_tail, head_tail_g
! 'This_side', the side the move comes from; 'that_side', the
! side the move goes to.
logical :: if_tangled
integer :: i, j, i0, j0, i1, j1, ib, jb, cv_nmb, cv_nmb_g
type(adss_info) :: new_ads, side_address
type(geo_info) :: gnb_orig
logical :: if_partner_delete

type(adss_info) :: next_stack
   
! type(geo_info) :: ggg_cell
   
!type(geo_info) :: g_g
   
! ggg_cell=temp%previous%g_cell
   
! type(auxi_crit_cell), pointer :: temp_xxx, temp_yyy
  
! nullify(temp_xxx); nullify(temp_yyy)
  
call move(temp,g_cell,p_cell,u_new)
! Produce the geometrical and physical information in the moved
! auxiliary critical cell.
i0=temp%g_cell%x_idx; j0=temp%g_cell%y_idx
i1=g_cell%x_idx; j1=g_cell%y_idx
! The indexes in the old and new auxiliary critical cells.
cross_o=-1000; cross_n=-1000
   
! The following sets the crossed corner points 'cross_o' and
! 'cross_n'.
select case(temp%g_cell%g_type)
 case('xx','yy')
  if(i1.gt.i0) then; cross_o(1)=2; cross_o(2)=3; end if
  if(i1.lt.i0) then; cross_o(1)=1; cross_o(2)=4; end if
  if(j1.gt.j0) then; cross_o(1)=3; cross_o(2)=4; end if
  if(j1.lt.j0) then; cross_o(1)=1; cross_o(2)=2; end if
  do j=1,2
   cross_n(j)=point_shift(cross_o(j),temp%g_cell%g_type)
  end do
 case('xy')
  if(i1.gt.i0) then
   if(j1.gt.j0) then
    cross_o(1)=3
   else
    cross_o(1)=2
   end if
  else
   if(j1.gt.j0) then
    cross_o(1)=4
   else
    cross_o(1)=1
   end if
  end if
  cross_n(1)=point_shift(cross_o(1),'xy ')
end select

! The following sets the sides.
if(move_dir(1).eq.'l_to_r') then
 this_side='left  '; that_side='right '
else
 this_side='right '; that_side='left  '
end if

! Update the spatial neighbor relation and the grid map as the 
! auxiliary critical cell leaves the old grid cell.
if(this_side.eq.'left  '.and.temp%l_stk%cv_nmb.le.0.or. &
   this_side.eq.'right '.and.temp%r_stk%cv_nmb.le.0) then
! The case that there is no other auxiliary critical cell stacked
! on 'this_side'.
 call clean(ggd_cell(i0,j0))
 uu(i0,j0)=u_new
else
! The case that there are other auxiliary critical cells stacked 
! on 'this_side'.
 if(this_side.eq.'left  ') then
  call visit(temp%l_stk,temp_g)
 else
  call visit(temp%r_stk,temp_g)
 end if
 do j=1,4
  if(ggd_cell(i0,j0)%ccpt(j)%address.eq.temp%address) then
   ggd_cell(i0,j0)%ccpt(j)%address=temp_g%address
   if(temp_g%l_stk.eq.temp%address) then
    ggd_cell(i0,j0)%ccpt(j)%side='left  '
   else
    if(temp_g%r_stk.eq.temp%address) then
     ggd_cell(i0,j0)%ccpt(j)%side='right '
    else
     call error_message
    end if
   end if
  end if
 end do
     		 	  	  	 	 	    
 if(temp_g%l_stk.eq.temp%address) then
  call clean_up_address(temp_g%l_stk)
 else
  call clean_up_address(temp_g%r_stk)
 end if

 if(this_side.eq.'left  ') then
  call clean_up_address(temp%l_stk)
 else
  call clean_up_address(temp%r_stk)
 end if
    
end if
   
! temp_xxx=>temp%previous 
! temp_yyy=>temp%next
   
if(temp%partner.ne.'previous') then
 call rm_tg(temp,g_cell,p_cell,this_side,that_side,if_tangled)
 if(if_tangled) then
  if_partner_delete=.false.
  if(temp%partner.eq.'next    ') then
   if(temp%next%g_cell%x_idx.eq.g_cell%x_idx.and.temp%next% &
      g_cell%y_idx.eq.g_cell%y_idx) then
    if_partner_delete=.true.
   else
    call error_message
   end if
  end if
  
  if(g_cell%g_type.ne.'xy ') then
   ib=temp%previous%g_cell%x_idx; jb=temp%previous%g_cell%y_idx
   if(i1.ne.ib.or.j1.ne.jb) then
    if(this_side.eq.'left  '.and.temp%l_stk.eq.adss_info(0,0)) then
     call clean(ggd_cell(i1,j1))
    end if
    if(this_side.eq.'right '.and.temp%r_stk.eq.adss_info(0,0)) then
     call clean(ggd_cell(i1,j1))
    end if
   end if
  end if
   
  call deletee(temp,'before')
  
  if(if_partner_delete) then
   
   if(temp%l_stk.eq.temp%next%address) then
    call clean_up_address(temp%l_stk)
   else
    if(temp%r_stk.eq.temp%next%address) then
     call clean_up_address(temp%r_stk)
    end if
   end if
    
   select case(move_dir(1))                                ! Modification Begin.
       
    case('l_to_r')  
          
	 if(temp%next%l_stk%cv_nmb.ne.0) then                   
      call visit(temp%next%l_stk,temp_gg)    
      temp%l_stk=temp_gg%address
      if(temp_gg%l_stk.eq.temp%next%address) then           
       temp_gg%l_stk=temp%address                           
      else
       if(temp_gg%r_stk.eq.temp%next%address) then
        temp_gg%r_stk=temp%address
       else
        call error_message
       end if
      end if
     end if    
	 
    case('r_to_l')  
          
     if(temp%next%r_stk%cv_nmb.ne.0) then                   
      call visit(temp%next%r_stk,temp_gg)    
      temp%r_stk=temp_gg%address
      if(temp_gg%l_stk.eq.temp%next%address) then           
       temp_gg%l_stk=temp%address                           
      else
       if(temp_gg%r_stk.eq.temp%next%address) then
        temp_gg%r_stk=temp%address
       else
        call error_message
       end if
      end if
     end if   	
     
    case default; call error_message
     
   end select                                             ! Modification End.	   		 	                                               
   
   temp=>temp%next
   call deletee(temp,'before') 
  end if
  return
 end if
end if
   
! Update the spatial neighbor relation and the grid map as the
! moved auxiliary critical cell arrives in the new grid cell.
if(temp%partner.eq.'single') then
! The case of no partner.
 if(ggd_cell(i1,j1)%region.eq.'smth') then
! The case that there is no other auxiliary critical cell to be
! stacked with the moved one in the new grid cell.
  uu(i1,j1)=error_data
  ggd_cell(i1,j1)%region='crit'
  do j=1,4
   ggd_cell(i1,j1)%ccpt(j)%address=temp%address
   if(j.eq.cross_n(1).or.j.eq.cross_n(2)) then
    ggd_cell(i1,j1)%ccpt(j)%side=this_side
   else
    ggd_cell(i1,j1)%ccpt(j)%side=that_side
   end if
  end do
 else
! The case there are other auxiliary critical cell to be stacked
! with the moved one in the new grid cell.
  call visit(ggd_cell(i1,j1)%ccpt(cross_n(1))%address,temp_g)

!  g_g=temp_g%g_cell
  if(temp%g_cell%g_type.ne.'xy ') then
   do j=1,2
    ggd_cell(i1,j1)%ccpt(cross_n(j))%address=temp%address
    ggd_cell(i1,j1)%ccpt(cross_n(j))%side=this_side
   end do
   if(that_side.eq.'left  ') temp%l_stk=temp_g%address
   if(that_side.eq.'right ') temp%r_stk=temp_g%address
  else
   temp%tangled='yes'
   ggd_cell(i1,j1)%ccpt(cross_n(1))%address=temp%address
   ggd_cell(i1,j1)%ccpt(cross_n(1))%side=that_side
   if(this_side.eq.'left  ') temp%l_stk=temp_g%address
   if(this_side.eq.'right ') temp%r_stk=temp_g%address
  end if

  head_tail='middle'; head_tail_g='middle'
  cv_nmb=temp%address%cv_nmb; cv_nmb_g=temp_g%address%cv_nmb
  if(.not.associated(temp%previous).and. &
     acvv(cv_nmb)%begin_end.gt.0) then
   if(.not.associated(temp_g%previous).and.acvv(cv_nmb)% &
      begin_end.eq.acvv(cv_nmb_g)%begin_end) then	  
    head_tail='begin '; head_tail_g='begin '
   end if
   if(.not.associated(temp_g%next).and.acvv(cv_nmb)% &
      begin_end.eq.acvv(cv_nmb_g)%end_end) then	  
    head_tail='begin '; head_tail_g='end   '
   end if
  end if
  if(.not.associated(temp%next).and.acvv(cv_nmb)%end_end.gt.0) then
   if(.not.associated(temp_g%previous).and.acvv(cv_nmb)% &
      end_end.eq.acvv(cv_nmb_g)%begin_end) then	  
    head_tail='end   '; head_tail_g='begin '
   end if
   if(.not.associated(temp_g%next).and.acvv(cv_nmb)% &
      end_end.eq.acvv(cv_nmb_g)%end_end) then	  
    head_tail='end   '; head_tail_g='end   '
   end if
  end if

  call clean_up_g_cell(gnb_orig)
  select case(head_tail)
   case('begin ','end   ')
    call side_neighbor(temp,temp_g,head_tail,head_tail_g)
   case('middle')
    call find_original(temp_g,gnb_orig)	   
    if(gnb_orig%point(cross_n(1)).eq.'left  ') then
     temp_g%l_stk=temp%address
    else
     temp_g%r_stk=temp%address
    end if
  end select
   
 end if          
 temp%g_cell=g_cell; temp%p_cell=p_cell; temp%updated='yes'
else
! The case of partners.
 nullify(temp_partner)
 select case(temp%partner)
  case('previous')
   temp%previous%partner='single  '; new_ads=temp%previous%address
   temp_partner=>temp%previous
  case('next    ')
   temp%next%partner='single  '; new_ads=temp%next%address
   temp_partner=>temp%next
 end select
 call deletee(temp,'before')
  
!  ggg_cell=temp%g_cell

 do i=1,2
  if(ggd_cell(i1,j1)%ccpt(cross_n(i))%address.eq.new_ads.and.ggd_cell(i1,j1)% &
     ccpt(cross_n(i))%side.eq.that_side) then
   ggd_cell(i1,j1)%ccpt(cross_n(i))%side=this_side
   ggd_cell(i1,j1)%ccpt(cross_n(i))%address=new_ads
  end if
  select case(that_side)
   case('left  '); side_address=temp_partner%l_stk
   case('right '); side_address=temp_partner%r_stk
   case default; call error_message
  end select    
  if(ggd_cell(i1,j1)%ccpt(cross_n(i))%address.eq.side_address) then
   ggd_cell(i1,j1)%ccpt(cross_n(i))%side=this_side
   ggd_cell(i1,j1)%ccpt(cross_n(i))%address=new_ads
  end if
 end do

end if    

end subroutine move_update


subroutine split_update(num)
! This subroutine updates the auxiliary critical cell that should
! split and the spatial neighbor relation and grid-map around the 
! cell.

implicit none

integer, intent(in) :: num
! The number of the discontinity position causing the split.

type(geo_info), dimension(2) :: gg_cell, gg_new, g_cell
type(phy_info), dimension(2) :: pp_cell, pp_new
! Geometrical and physical information in the original and
! split auxiliary critical cells. If there is partnership, 
! two pairs of them will be involved.
integer :: cross_o, cross_n
! The corner point the split will cross, 'cross_o', viewed in
! the old cell and 'cross_n', viewed in the new cell.
type(auxi_crit_cell), pointer :: temp_g, temp_new
! Pointer to the neighboring auxiliary critical cell and the new
! critical cell generated in the split.
character*6 :: this_side, that_side, head_tail, head_tail_g
! 'This_side', the side the move comes from; 'that_side', the
! side the move goes to.

type(geo_info) :: gnb_orig, gnn_cell, gg_c
integer, dimension(2) :: numb
character*6, dimension(2,2) :: mv_dir
integer :: involved, i, j, i0, j0, i1, j1, cv_nmb, cv_nmb_g
character*3 :: direction, cross_direction
logical :: if_tangled, reverse
type(state) :: difference
integer :: corner_crossed, point, single, ii
type(adss_info) :: adss_partner, ttt
type(auxi_crit_cell), pointer :: temp_neighb

integer :: ssingle

! type(adss_info) :: address_tran
! type(geo_info) :: ggg

nullify(temp_g); reverse=.false.
g_cell(1)=temp%g_cell
gg_cell(1)=temp%g_cell; pp_cell=temp%p_cell
involved=1; numb(1)=num; numb(2)=-1000
mv_dir(1,:)=move_dir; mv_dir(2,:)=' '
call split(temp,gg_cell(1),gg_new(1),pp_cell(1),pp_new(1),numb(1))
! Geometrically and physically split the auxiliary critical cell
! first.

if(gg_cell(1)%x_idx.ne.gg_new(1)%x_idx) direction='xx '
if(gg_cell(1)%y_idx.ne.gg_new(1)%y_idx) direction='yy '
! Determine the direction of split.

i0=gg_cell(1)%x_idx; j0=gg_cell(1)%y_idx
if(temp%partner.eq.'next    ') then
 temp_g=>temp%next
 i1=temp_g%g_cell%x_idx; j1=temp_g%g_cell%y_idx

 if(gg_new(1)%x_idx.eq.i1.and.gg_new(1)%y_idx.eq.j1) then
  if(num.eq.1) reverse=.true.
 end if
 if(gg_new(1)%x_idx.eq.i1.and.gg_new(1)%y_idx.eq.j1) then
! If there is partner relation, the partner should also be split.
  
  involved=2; numb(2)=3-num
  if(mv_dir(1,numb(1)).eq.'l_to_r') mv_dir(2,numb(2))='r_to_l'
  if(mv_dir(1,numb(1)).eq.'r_to_l') mv_dir(2,numb(2))='l_to_r'
  move_dir=mv_dir(2,:)
  g_cell(2)=temp_g%g_cell
  gg_cell(2)=temp_g%g_cell; pp_cell(2)=temp_g%p_cell
  call split(temp_g,gg_cell(2),gg_new(2),pp_cell(2),pp_new(2),numb(2))
! Adjust the geometrical and physical information for the pair of
! partnered auxiliary critical cells.
  pp_cell(1)%or_state=0.5d0*(pp_cell(1)%or_state+pp_new(2)%or_state)
  pp_cell(2)%or_state=0.5d0*(pp_cell(2)%or_state+pp_new(1)%or_state)
  if(.not.reverse) then
   gg_cell(1)%dis_posi(2)=0.5d0* &
   (gg_cell(1)%dis_posi(2)+gg_cell(2)%dis_posi(1)) 
   gg_cell(2)%dis_posi(1)=gg_cell(1)%dis_posi(2)
  else
   gg_cell(1)%dis_posi(1)=0.5d0* &
   (gg_cell(1)%dis_posi(1)+gg_cell(2)%dis_posi(2)) 
   gg_cell(2)%dis_posi(2)=gg_cell(1)%dis_posi(1)
  end if 
 end if
end if

do i=1,involved 
 i0=gg_cell(i)%x_idx; j0=gg_cell(i)%y_idx
 i1=gg_new(i)%x_idx; j1=gg_new(i)%y_idx
 cross_o=-1000; cross_n=-1000
! Set the corssed corner point.
 if(i1+j1.gt.i0+j0) then
  select case(g_cell(i)%edge(numb(i)))
   case(1,2); cross_o=cycle_4(g_cell(i)%edge(numb(i))+1)
   case(3,4); cross_o=g_cell(i)%edge(numb(i))
  end select
 else
  select case(g_cell(i)%edge(numb(i)))
   case(1,2); cross_o=g_cell(i)%edge(numb(i))
   case(3,4); cross_o=cycle_4(g_cell(i)%edge(numb(i))+1)
  end select
 end if
 cross_n=point_shift(cross_o,direction)
! Set the sides.
 if(mv_dir(i,numb(i)).eq.'l_to_r') then
  this_side='left  '; that_side='right '
 else
  this_side='right '; that_side='left  '
 end if

!  ggg=temp%g_cell

 call clean_up_address(adss_partner)
 select case(this_side)
  case('left  ')
  if(i.eq.1) then
   adss_partner=temp%r_stk
  else 
   adss_partner=temp_g%r_stk
  end if
  case('right ')
   if(i.eq.1) then
    adss_partner=temp%l_stk
   else
    adss_partner=temp_g%l_stk
   end if
  case default; call error_message
 end select

! Update of the grid map caused by departure.
 if(i.eq.1) then
  if(ggd_cell(i0,j0)%ccpt(cross_o)%address.eq.temp%address) then
   ggd_cell(i0,j0)%ccpt(cross_o)%side=this_side
  else
   if(ggd_cell(i0,j0)%ccpt(cross_o)%address.eq.adss_partner) then
    ggd_cell(i0,j0)%ccpt(cross_o)%address=temp%address
    ggd_cell(i0,j0)%ccpt(cross_o)%side=this_side
   end if
  end if
 else
  if(ggd_cell(i0,j0)%ccpt(cross_o)%address.eq.temp_g%address) then
   ggd_cell(i0,j0)%ccpt(cross_o)%side=this_side
  else
   if(ggd_cell(i0,j0)%ccpt(cross_o)%address.eq.adss_partner) then
    ggd_cell(i0,j0)%ccpt(cross_o)%address=temp_g%address
    ggd_cell(i0,j0)%ccpt(cross_o)%side=this_side
   end if
  end if
 end if 
end do

! Update of the spatial neighbor relation and grid map caused
! by arrival.
if(involved.eq.1) then
! No partner relation or partnered with te next critical cell.
 temp%g_cell=gg_cell(1); temp%p_cell=pp_cell(1); temp%updated='yes'
 allocate(temp_new)
 nullify(temp_new%p_nxt_dr); nullify(temp_new%n_nxt_dr)
 call clean_up_address(temp_new%l_stk); call clean_up_address(temp_new%r_stk)
 temp_new%g_cell=gg_new(1); temp_new%p_cell=pp_new(1)
 temp_new%address%cv_nmb=temp%address%cv_nmb
 temp_new%updated='yes'; temp_new%partner='single  ' 
 temp_new%stencil=error_data; temp_new%stencil_disp=error_index
! Generate the pointer pointing to the split auxiliary critical 
! cell.

  temp%p_cell%or_state%gamma=error_data
  temp_new%p_cell%or_state%gamma=error_data

 if(numb(1).eq.1) then
  call insert(temp,temp_new,'before')
 else
  call insert(temp,temp_new,'after ')
 end if

 if(numb(1).eq.1) then
  call rm_tg(temp_new,gg_new(1),pp_new(1),this_side,that_side,if_tangled)
  if(if_tangled) then
   call deletee(temp_new,'before'); return
  end if
 end if  

 if(ggd_cell(i1,j1)%region.eq.'smth') then
! No other auxiliary critical cell to be stacked with the split
! one in the new grid cell.
  uu(i1,j1)=error_data
  ggd_cell(i1,j1)%region='crit'
  do j=1,4
   ggd_cell(i1,j1)%ccpt(j)%address=temp_new%address
   if(j.eq.cross_n) then
    ggd_cell(i1,j1)%ccpt(j)%side=this_side
   else
    ggd_cell(i1,j1)%ccpt(j)%side=that_side
   end if
  end do
 else
! There are other auxiliary critical cells to be stacked with
! the split one in the new grid cell.
  call visit(ggd_cell(i1,j1)%ccpt(cross_n)%address,temp_g)
  ggd_cell(i1,j1)%ccpt(cross_n)%address=temp_new%address
  ggd_cell(i1,j1)%ccpt(cross_n)%side=this_side
  if(that_side.eq.'left  ') temp_new%l_stk=temp_g%address
  if(that_side.eq.'right ') temp_new%r_stk=temp_g%address

  head_tail='middle'; head_tail_g='middle'
  cv_nmb=temp%address%cv_nmb; cv_nmb_g=temp_g%address%cv_nmb
  if(.not.associated(temp_new%previous).and. &
     acvv(cv_nmb)%begin_end.gt.0.and.num.eq.1) then
   if(.not.associated(temp_g%previous).and.acvv(cv_nmb)% &
      begin_end.eq.acvv(cv_nmb_g)%begin_end) then
    head_tail='begin '; head_tail_g='begin '
   end if
   if(.not.associated(temp_g%next).and.acvv(cv_nmb)% &
      begin_end.eq.acvv(cv_nmb_g)%end_end) then
    head_tail='begin '; head_tail_g='end   '
   end if
  end if
  if(.not.associated(temp_new%next).and. &
     acvv(cv_nmb)%end_end.gt.0.and.num.eq.2) then
   if(.not.associated(temp_g%previous).and.acvv(cv_nmb)% &
      end_end.eq.acvv(cv_nmb_g)%begin_end) then
    head_tail='end   '; head_tail_g='begin '
   end if
   if(.not.associated(temp_g%next).and.acvv(cv_nmb)% &
      end_end.eq.acvv(cv_nmb_g)%end_end) then
    head_tail='end   '; head_tail_g='end   '
   end if
  end if

! For the two ends of the discontinuity curve, special care
! should be taken in updating the spatial relation.
  call clean_up_g_cell(gnb_orig)
  select case(head_tail)
   case('begin ','end   ')
    call side_neighbor(temp_new,temp_g,head_tail,head_tail_g)
   case('middle')
    call find_original(temp_g,gnb_orig)	   
    if(gnb_orig%point(cross_n).eq.'left  ') then
     temp_g%l_stk=temp_new%address
    else
     temp_g%r_stk=temp_new%address
    end if
  end select

 end if
 
 if(num.eq.1.and.temp%partner.eq.'next  ') then
!  if(gg_cell(1)%g_type.eq.'xy ') then
   gg_c=temp%g_cell
   gnn_cell=temp%next%g_cell
   if(dabs(gnn_cell%dis_posi(2)).lt.0.5d0) call error_message
! The following determine the corner point the discontinity is 
! going to cross.
   corner_crossed=error_index
   select case(gnn_cell%edge(2))
    case(1)
     if(gg_c%x_idx.gt.gnn_cell%x_idx) corner_crossed=2
     if(gg_c%x_idx.lt.gnn_cell%x_idx) corner_crossed=1
    case(2)
     if(gg_c%y_idx.gt.gnn_cell%y_idx) corner_crossed=3
     if(gg_c%y_idx.lt.gnn_cell%y_idx) corner_crossed=2
    case(3)
     if(gg_c%x_idx.gt.gnn_cell%x_idx) corner_crossed=3
     if(gg_c%x_idx.lt.gnn_cell%x_idx) corner_crossed=4
    case(4)
     if(gg_c%y_idx.gt.gnn_cell%y_idx) corner_crossed=4
     if(gg_c%y_idx.lt.gnn_cell%y_idx) corner_crossed=1
   end select
   if(corner_crossed.eq.error_index) call error_message

 ! The following determine the corss direction. 
   cross_direction='   '
   if(gg_c%x_idx.ne.gnn_cell%x_idx) cross_direction='xx '
   if(gg_c%y_idx.ne.gnn_cell%y_idx) cross_direction='yy '
   if(cross_direction.eq.'   ') call error_message

   point=point_shift(corner_crossed,cross_direction)
   	     
! The following updates the grid-map at this point.
   ggd_cell(gg_c%x_idx,gg_c%y_idx)%ccpt(point)%address= &
   temp%address
   ggd_cell(gg_c%x_idx,gg_c%y_idx)%ccpt(point)%side= &
   gg_c%point(point)
     
   temp=>temp%next
   call remove(temp,'before',gnn_cell%point(corner_crossed), &
               difference)
!  end if 
 end if

 if(numb(1).ne.1) temp=>temp%next 
else
! Partner relation.
 if(.not.reverse) then
  if(numb(1).eq.1) then
   temp%g_cell=gg_cell(2); temp%p_cell=pp_cell(2)
   temp_g%g_cell=gg_cell(1); temp_g%p_cell=pp_cell(1)
  else
   temp%g_cell=gg_cell(1); temp%p_cell=pp_cell(1)
   temp_g%g_cell=gg_cell(2); temp_g%p_cell=pp_cell(2)
  endif
 else

   print*, '   '
   print*, '   '
   print*,'Reverse happens at i=',temp%g_cell%x_idx, 'and j=', temp%g_cell%y_idx
   print*, '   '
   print*, '   '

  temp%g_cell=gg_cell(2); temp%p_cell=pp_cell(2)
  temp_g%g_cell=gg_cell(1); temp_g%p_cell=pp_cell(1)

  i0=gg_cell(2)%x_idx; j0=gg_cell(2)%y_idx
  do ii=1, 4
   if(ggd_cell(i0,j0)%ccpt(ii)%address.eq.temp_g%address) then
    ggd_cell(i0,j0)%ccpt(ii)%address=temp%address
    ggd_cell(i0,j0)%ccpt(ii)%side=gg_cell(2)%point(ii)
   end if
  end do   
  
!%%%%%%%%%%%%%%%%%% correction %%%%%%%%%

  call find_single(gg_cell(2),ssingle)
  ggd_cell(i0,j0)%ccpt(ssingle)%address=temp%address
  ggd_cell(i0,j0)%ccpt(ssingle)%side=gg_cell(1)%point(ssingle) 

!%%%%%%%%%%%%%%%%%% correction %%%%%%%%%
    
  i0=gg_cell(1)%x_idx; j0=gg_cell(1)%y_idx
  do ii=1, 4
   if(ggd_cell(i0,j0)%ccpt(ii)%address.eq.temp%address) then
    ggd_cell(i0,j0)%ccpt(ii)%address=temp_g%address
	ggd_cell(i0,j0)%ccpt(ii)%side=gg_cell(1)%point(ii)
   end if
  end do 
  
!%%%%%%%%%%%%%%%%%% correction %%%%%%%%%

  call find_single(gg_cell(1),ssingle)
  ggd_cell(i0,j0)%ccpt(ssingle)%address=temp_g%address
  ggd_cell(i0,j0)%ccpt(ssingle)%side=gg_cell(1)%point(ssingle) 

!%%%%%%%%%%%%%%%%%% correction %%%%%%%%%
  

! If reserve situation happens, the corresponding neighboring relations 
! should be updated.

! First, update both the left neighboring relations.
  ttt=temp%l_stk; temp%l_stk=temp_g%l_stk; temp_g%l_stk=ttt
  
  nullify(temp_neighb)
  if(temp_g%l_stk%cv_nmb.ne.0) then
   call visit(temp_g%l_stk,temp_neighb)
   if(temp_neighb%l_stk.eq.temp%address) then
    temp_neighb%l_stk=temp_g%address
   else
    if(temp_neighb%r_stk.eq.temp%address) then
     temp_neighb%r_stk=temp_g%address
    else
	 call error_message
    end if
   end if		  
  end if
    
  nullify(temp_neighb)
  if(temp%l_stk%cv_nmb.ne.0) then
   call visit(temp%l_stk,temp_neighb)
   if(temp_neighb%l_stk.eq.temp_g%address) then
    temp_neighb%l_stk=temp%address
   else
    if(temp_neighb%r_stk.eq.temp_g%address) then
     temp_neighb%r_stk=temp%address
    else
	 call error_message
    end if
   end if		  
  end if 

! Second, update both the right neighboring relations.
  ttt=temp%r_stk; temp%r_stk=temp_g%r_stk; temp_g%r_stk=ttt

  nullify(temp_neighb)
  if(temp_g%r_stk%cv_nmb.ne.0) then
   call visit(temp_g%r_stk,temp_neighb)
   if(temp_neighb%l_stk.eq.temp%address) then
    temp_neighb%l_stk=temp_g%address
   else
    if(temp_neighb%r_stk.eq.temp%address) then
     temp_neighb%r_stk=temp_g%address
    else
	 call error_message
    end if
   end if		  
  end if
    
  nullify(temp_neighb)
  if(temp%r_stk%cv_nmb.ne.0) then
   call visit(temp%r_stk,temp_neighb)
   if(temp_neighb%l_stk.eq.temp_g%address) then
    temp_neighb%l_stk=temp%address
   else
    if(temp_neighb%r_stk.eq.temp_g%address) then
     temp_neighb%r_stk=temp%address
    else
	 call error_message
    end if
   end if		  
  end if 

!   call find_single(gg_cell(ii),single)
!   i0=gg_cell(ii)%x_idx; j0=gg_cell(ii)%y_idx
!   ggd_cell(i0,j0)%ccpt(single)%side=side_shift(ggd_cell(i0,j0) &
!                                     %ccpt(single)%side)
 end if

 temp%partner='single  '; temp_g%partner='single  '
 temp_g%updated='yes'

end if    
temp%updated='yes'

move_dir=mv_dir(1,:); move_dir(num)=' '

end subroutine split_update


subroutine fix_heads(acv)
! This subroutine fixes the tangle of the begin-end and end-end.

implicit none
type(auxi_discv), intent(inout) :: acv
! The single auxiliary discontinuity curve to be updated.

type(geo_info) :: g_cell_p
type(phy_info) :: p_cell_p
integer :: i0, j0, i1, j1
character*6 :: this_side, that_side
logical :: if_tangled

type(auxi_crit_cell), pointer :: temp_p
nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
end if

g_cell=temp%g_cell; p_cell=temp%p_cell
i0=g_cell%x_idx; j0=g_cell%y_idx
temp_p=>temp%previous
g_cell_p=temp_p%g_cell; p_cell_p=temp_p%p_cell
i1=g_cell_p%x_idx; j1=g_cell_p%y_idx

if(i1.eq.i0.and.j1.eq.j0) then
 if(temp_p%l_stk.eq.temp%address) then
  this_side='right '; that_side='left  '
 else
  if(temp_p%r_stk.eq.temp%address) then
   this_side='left  '; that_side='right '
  else
   call error_message
  end if
 end if
 temp%l_stk=adss_info(0,0); temp%r_stk=adss_info(0,0)
 temp_p%l_stk=adss_info(0,0); temp_p%r_stk=adss_info(0,0)
  
 call rm_tg(temp,g_cell,p_cell,this_side,that_side,if_tangled)
 if(if_tangled) then
  call deletee(temp,'before'); return
 end if

end if

end subroutine fix_heads


subroutine fix_remaining_tangled(acv,possible_remaining)
! This subroutine fixes the tangle of the begin-end and end-end.

implicit none
type(auxi_discv), intent(inout) :: acv
! The single auxiliary discontinuity curve to be updated.
logical, intent(inout) :: possible_remaining

type(geo_info) :: g_cell_p
type(phy_info) :: p_cell_p
integer :: i0, j0, i1, j1, i
character*6 :: side_between, side_apart
type(auxi_crit_cell), pointer :: temp_p, temp_g, temp_g1, boundary_end
type(state) :: uuu, or_uu, difference

possible_remaining=.false.

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
  
 g_cell=temp%g_cell; p_cell=temp%p_cell
 i0=g_cell%x_idx; j0=g_cell%y_idx
 temp_p=>temp%previous
 g_cell_p=temp_p%g_cell; p_cell_p=temp_p%p_cell
 i1=g_cell_p%x_idx; j1=g_cell_p%y_idx
  
 if(i1.eq.i0.and.j1.eq.j0) then
! There is a remaining tangled in the cell and the two tangled discontinuity
! cells are facing to each other.

! We need to determine the side between and the side on the two sides of 
! the two tangled discontinuity cells.
  if(temp_p%l_stk.eq.temp%address) then
   if(temp%l_stk.eq.temp_p%address) then
    side_apart='right '; side_between='left  '
   else
    call error_message
   end if
  else
   if(temp_p%r_stk.eq.temp%address) then
    if(temp%r_stk.eq.temp_p%address) then
     side_apart='left  '; side_between='right '
    else
     call error_message
    end if
   else
    call error_message
   end if
  end if
  
  if(g_cell%edge(2).ne.g_cell_p%edge(1)) then
! In this case, the two discontinuity cells should be merged.
   
! Firstly, update the previous discontinuity cell, which will become the
! merged discontinuity cell and thus remain.
   p_cell_p%l_state=0.5d0*(p_cell%l_state+p_cell_p%l_state)
   p_cell_p%r_state=0.5d0*(p_cell%r_state+p_cell_p%r_state)
   p_cell_p%or_state=p_cell%or_state+p_cell_p%or_state
   select case(side_between)
    case('left  '); p_cell_p%or_state=p_cell_p%or_state-p_cell_p%l_state
    case('right '); p_cell_p%or_state=p_cell_p%or_state-p_cell_p%r_state
   end select
    
! Modification 2016-3-10, bounds for `or_state' when 'remaining tangled' is treated.
  do i=1, 4
   p_cell_p%or_state%value(i)=between(p_cell_p%or_state%value(i),p_cell_p%l_state%value(i),p_cell_p%r_state%value(i))   
  end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  temp_p%p_cell=p_cell_p
  g_cell_p%edge(2)=g_cell%edge(2)
  g_cell_p%dis_posi(2)=g_cell%dis_posi(2)
  call find_type_from_edges(g_cell_p)
  call find_points_from_edges(g_cell_p)       
  temp_p%g_cell=g_cell_p
   
! Secondly, updated the grid_map.     
  do i=1,4
   if(ggd_cell(i1,j1)%ccpt(i)%address.eq.temp%address) then
    ggd_cell(i1,j1)%ccpt(i)%address=temp_p%address
    ggd_cell(i1,j1)%ccpt(i)%side=g_cell_p%point(i)
   end if 
  end do 
   
! Thirdly, update the neighboring relation.
   nullify(temp_g)
   select case(side_apart)
    case('left  ')
     if(temp%l_stk%cv_nmb.ne.0.and.temp_p%l_stk%cv_nmb.ne.0) then
      call error_message
     end if
     if(temp%l_stk%cv_nmb.ne.0) then
      temp_p%r_stk=temp%l_stk
      call visit(temp%l_stk,temp_g)
     else
      temp_p%r_stk=adss_info(0,0)	  
     end if
    case('right ')
     if(temp%r_stk%cv_nmb.ne.0.and.temp_p%r_stk%cv_nmb.ne.0) then
      call error_message
     end if
     if(temp%r_stk%cv_nmb.ne.0)  then
      temp_p%l_stk=temp%r_stk
      call visit(temp%r_stk,temp_g)
     else
	  temp_p%l_stk=adss_info(0,0)
     end if 
   end select
   
   if(associated(temp_g)) then
    if(temp_g%l_stk.eq.temp%address) then
     temp_g%l_stk=temp_p%address
    else
     if(temp_g%r_stk.eq.temp%address) then
      temp_g%r_stk=temp_p%address
     else  
      call error_message
     end if
    end if
   end if

 ! Finally, delete the current discontinuity cell.
   call deletee(temp,'before')
   
!    call check_list_ac(1,'down  ')
  else
! In this case, the two discontinuity cells should be deleted.
   
! Firstly, update the solution in the cell, which is now in smooth region.
   select case(side_apart)
    case('left  ')
     uuu=0.5d0*(p_cell%l_state+p_cell_p%l_state)
     or_uu=p_cell%or_state+p_cell_p%or_state-p_cell%r_state
    case('right ')
     uuu=0.5d0*(p_cell%r_state+p_cell_p%r_state)
     or_uu=p_cell%or_state+p_cell_p%or_state-p_cell%l_state
   end select
   difference=or_uu-uuu   
   temp_p%previous%p_cell%or_state=temp_p%previous%p_cell%or_state &
                                   +difference
   temp%next%p_cell%or_state=temp%next%p_cell%or_state+difference   
   
! Secondly, update the grid_map and the neighboring relation.
   nullify(temp_g); nullify(temp_g1)
   select case(side_apart)
    case('left  ')
    if(temp_p%l_stk%cv_nmb.eq.0.and.temp%l_stk%cv_nmb.eq.0) then
      call clean(ggd_cell(i1,j1))
      uu(i1,j1)=uuu
      goto 10
     else
      if(temp%l_stk%cv_nmb.ne.0) call visit(temp%l_stk,temp_g)
      if(temp_p%l_stk%cv_nmb.ne.0) call visit(temp_p%l_stk,temp_g1)
     end if
    case('right ')
     if(temp_p%r_stk%cv_nmb.eq.0.and.temp%r_stk%cv_nmb.eq.0) then
      call clean(ggd_cell(i1,j1))
      uu(i1,j1)=uuu
      goto 10
     else
      if(temp%r_stk%cv_nmb.ne.0) call visit(temp%r_stk,temp_g)
      if(temp_p%r_stk%cv_nmb.ne.0) call visit(temp_p%r_stk,temp_g1)
     end if	
   end select	   	       
   
   if(associated(temp_g)) then
    if(temp_g%l_stk.eq.temp%address) then 
     do i=1,4 
      if(ggd_cell(i1,j1)%ccpt(i)%address.eq.temp%address.or. &
	     ggd_cell(i1,j1)%ccpt(i)%address.eq.temp_p%address) then
       ggd_cell(i1,j1)%ccpt(i)%address=temp_g%address
       ggd_cell(i1,j1)%ccpt(i)%side='left  '
      end if
     end do
     if(associated(temp_g1)) then	  
      temp_g%l_stk=temp_g1%address
     else
      temp_g%l_stk=adss_info(0,0)
     end if	  	  
    else
     if(temp_g%r_stk.eq.temp%address) then
      do i=1,4 
       if(ggd_cell(i1,j1)%ccpt(i)%address.eq.temp%address.or. &
	      ggd_cell(i1,j1)%ccpt(i)%address.eq.temp_p%address) then
        ggd_cell(i1,j1)%ccpt(i)%address=temp_g%address
        ggd_cell(i1,j1)%ccpt(i)%side='right '
       end if
      end do
      if(associated(temp_g1)) then	   
       temp_g%r_stk=temp_g1%address
      else
       temp_g%r_stk=adss_info(0,0)
      end if
     else
      call error_message
     end if
    end if
   end if

   if(associated(temp_g1)) then        
    if(temp_g1%l_stk.eq.temp_p%address) then 
     do i=1,4 
      if(ggd_cell(i1,j1)%ccpt(i)%address.eq.temp_p%address.or. &
	     ggd_cell(i1,j1)%ccpt(i)%address.eq.temp%address) then
       ggd_cell(i1,j1)%ccpt(i)%address=temp_g1%address
       ggd_cell(i1,j1)%ccpt(i)%side='left  '
      end if
     end do 
     if(associated(temp_g)) then
      temp_g1%l_stk=temp_g%address
     else
      temp_g1%l_stk=adss_info(0,0)
     end if	   
    else
     if(temp_g1%r_stk.eq.temp_p%address) then
      do i=1,4 
       if(ggd_cell(i1,j1)%ccpt(i)%address.eq.temp_p%address.or. &
	      ggd_cell(i1,j1)%ccpt(i)%address.eq.temp%address) then
        ggd_cell(i1,j1)%ccpt(i)%address=temp_g1%address
        ggd_cell(i1,j1)%ccpt(i)%side='right '
       end if
      end do 
      if(associated(temp_g)) then
       temp_g1%r_stk=temp_g%address
      else
       temp_g1%r_stk=adss_info(0,0)
      end if
     else
      call error_message
     end if
    end if
   end if
! Finally, delete the two discotinuity cells.
   
10 call deletee(temp_p,'before')
   call deletee(temp,'before')
   possible_remaining=.true.
    
  end if
     
 end if
   
 temp=>temp%next
 if(associated(temp)) then
 head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
 else
  exit
 end if
  
end do
  
end subroutine fix_remaining_tangled


subroutine head_shift(acv)
! This subroutine checks that the two heads of the curve are at a
! right place so that no complications will be caused in the reset.
! If the heads are found not at a right place, this subroutine will
! then shift them to a right place.

implicit none
type(auxi_discv), intent(inout) :: acv
! The single auxiliary discontinuity curve to be updated.

type(geo_info) :: g_cell_p
logical :: if_shift

nullify(temp)
if(associated(acv%begin%next)) then
 temp=>acv%begin%next
 head_mark=temp%address%idx; head_switch=.true.
end if

if_shift=.false.
g_cell=temp%g_cell
g_cell_p=temp%previous%g_cell

if(g_cell%g_type.eq.'xy '.or.g_cell%g_type.eq.'xy ') then
 if_shift=.true.
else
 if(dabs(g_cell%dis_posi(1)).gt.0.5d0.or.dabs( &
    g_cell%dis_posi(2)).gt.0.5d0) if_shift=.true.
end if

if(.not.if_shift) return

do while(associated(temp).and.head_switch) 

 if_shift=.false.
 g_cell=temp%g_cell
 g_cell_p=temp%previous%g_cell

 if(g_cell%g_type.eq.'xy '.or.g_cell%g_type.eq.'xy ') then
  if_shift=.true.
 else
  if(dabs(g_cell%dis_posi(1)).gt.0.5d0.or.dabs( &
     g_cell%dis_posi(2)).gt.0.5d0) if_shift=.true.
 end if

 if(.not.if_shift) then
  acv%begin%next=>temp
  acv%eend%previous=>temp%previous
  return
 end if

 temp=>temp%next
 if(associated(temp)) then
 head_switch=(temp%address%idx.ne.head_mark)
 else
  exit
 end if
  
end do

print*, 'There is no match in the head-shift!' 
call error_message

end subroutine head_shift


end module update_curves
