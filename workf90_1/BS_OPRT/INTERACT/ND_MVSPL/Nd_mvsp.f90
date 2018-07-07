module node_move_n_split
! This module deals with move and split of node cells

use merge_checks
! 'merge_ck.f90'

use node_production_n_update_1
! 'nd_proc1.f90'

use node_production_n_update_2
! 'nd_proc2.f90'

use connect_and_triple
! 'contrip.f90'

use corner_remove
! 'rm_cons.f90'

use disjoints_fix
! 'dis_jntn.f90'

public move_n_splits
private uu, cvv, ggd_cell, ndd, check1, check2, check3, &
        nd_nmb, nd_nmb_n, x_nd, y_nd, x_nd_new, y_nd_new, &
		xw_posi, yw_posi, if_merge, current_merged, next_merged, &
		previous_merged, tempp, tempp_n, head_mark, head_switch, &
		plugs_update,  merge_check, find_edgec, sum, &
        evaluating_node_state, find_node, locate_head, &
		find_node_number, find_head_cr, find_plug_posi


contains


subroutine move_n_splits

implicit none
 
! integer :: i0, j0

allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(cvv(curves_number))
allocate(ndd(ndn))

call give_sth(0)
call give_cv
call give_gcl
call give_nd

! call scanning(uu)

! call check_list_c(4,'up    ')

! call check_ring('counter ')
 
!  i0=ndd(1)%n_cell%x_idx; j0=ndd(1)%n_cell%y_idx

! call check_conservation(-4,0,4,5)

do nd_nmb=1, ndn
 if(ndd(nd_nmb)%status.eq.'awake ') then

!   call check_list_c(4,'up    ')
!   call check_ring('counter ')

  call move_n_split

!   call scanning(uu)
!   call check_ring('counter ')
!   call check_list_c(4,'down  ')

  call connect_n_triple

!   call check_list_c(4,'up    ')   

 end if
end do

! call check_conservation(-4,0,4,5)

! call check_list_c(4,'up    ')

! call check_ring(2,'counter ')

call waking_nodes
call waking_curves
call remove_corners
call fix_all_disjoints
call clean_node_cells

call get_sth(0)
call get_cv
call get_gcl
call get_nd

! call scanning(uu)

deallocate(uu)
deallocate(ggd_cell)
deallocate(cvv)
deallocate(ndd)

end subroutine move_n_splits


subroutine move_n_split
! This subroutine inplements the move and splitting of a node
! cell.

implicit none

integer :: i

type(cv_plug_info) :: plugg
logical :: physical_possibility
! type(state) :: xx

sum=ndd(nd_nmb)%n_cell%or_state
x_nd=ndd(nd_nmb)%n_cell%x_idx; y_nd=ndd(nd_nmb)%n_cell%y_idx
new_nodes=0; current_node_empty='no '
tempp=>ndd(nd_nmb)%cv_rg%begin%next
head_mark=tempp%address%idx; head_switch=.true.
do while(head_switch)

  plugg=tempp%plug
!  call check_ring(1,'counter ')
 if_merge='no '

 call if_physically_possible
 ! Check if the two discontinuities can physically interact with each other. For examples,
 ! in Euler system, contact discontinuities or slip lines can not catch a shock from behind, 
 ! or two shocks of different types can not collide each other from behind.

 if(physical_possibility) call check_mergence
! Check whether there is a mergence of discontinuity curves.

!  plugg=tempp_n%plug

!  call check_ring(1,'counter ')

 if(if_merge.ne.'no    ') then 
  new_nodes=new_nodes+1
  call plugs_update

!   call check_ring(1,'clock   ')
!   plugg=tempp%plug
!   xx=uu(-2,1)

  call connect_n_produce
   
!   call check_list_c(7,'up    ')
    
  if(current_node_empty.eq.'yes') exit
 end if
! Produce the new node cell caused by the mergence and update 
! the old node cell.

!  call check_ring(1,'clock   ')
!  plugg=tempp_n%plug

 if(associated(tempp)) then
  select case(if_merge)
   case('yes'); tempp=>tempp_n
   case('no '); tempp=>tempp%next
   case default; call error_message
  end select
!   plugg=tempp%plug  
  head_switch=(tempp%address%idx.ne.head_mark)
 else
  head_switch=.false.
 end if
end do

ndd(nd_nmb)%n_cell%or_state=sum
if(current_node_empty.eq.'no ') call renumber_plugs(nd_nmb)

! The following puts 'awake' empty discontinuity curve list back
! to 'yawn' status.
do i=1, cvn; call put_empty_yawn(i); end do

! call check_ring(1,'counter ')


contains


subroutine if_physically_possible

implicit none

type(curve_plug), pointer :: tempp_n
type(cv_plug_info) :: plugg_n
type(phy_info) :: p_cell_ph
type(critical_cell), pointer :: tempph
integer :: wave_nmb, wave_nmb_n, cvv_nmb, cvv_nmb_n

tempp_n=>tempp%next

plugg_n=tempp_n%plug

cvv_nmb=plugg%cv_nmb; cvv_nmb_n=plugg_n%cv_nmb
if(cvv(cvv_nmb)%status.ne.'awake '.or.cvv(cvv_nmb_n)%status.ne.'awake ') then
 physical_possibility=.false.
 return
end if  

select case(plugg%end_type)
 case('begin '); tempph=>cvv(plugg%cv_nmb)%begin%next
 case('end   '); tempph=>cvv(plugg%cv_nmb)%eend%previous
 case default; call error_message
end select
p_cell_ph=tempph%p_cell
select case(plugg%end_type)
 case('begin '); wave_nmb=p_cell_ph%wv_nb
 case('end   '); wave_nmb=wave_number+1-p_cell_ph%wv_nb
end select

select case(plugg_n%end_type)
 case('begin '); tempph=>cvv(plugg_n%cv_nmb)%begin%next
 case('end   '); tempph=>cvv(plugg_n%cv_nmb)%eend%previous
 case default; call error_message
end select
p_cell_ph=tempph%p_cell
select case(plugg_n%end_type)
 case('begin '); wave_nmb_n=p_cell_ph%wv_nb
 case('end   '); wave_nmb_n=wave_number+1-p_cell_ph%wv_nb
end select

if(wave_nmb.le.wave_nmb_n) then
 physical_possibility=.true.
else
 physical_possibility=.false.
end if

end subroutine if_physically_possible



end subroutine move_n_split


subroutine check_conservation(i0,j0,nd_nmb,nd_nmb_n)
! This subroutine is designed to check the conservation in the
! vicinity of the node cell(s).

implicit none

integer, intent(in) :: i0, j0, nd_nmb, nd_nmb_n

type(critical_cell), pointer :: temp, temp_n
type(adss_info) :: ads, ads_n
type(geo_info) :: g_cell, g_cell_n
type(phy_info) :: p_cell, p_cell_n
type(state), dimension(-1:1,-1:1) :: uuu
type(state) :: conservation
integer :: i, j

do i=-1,1; do j=-1,1
 uuu(i,j)=0.d0 
end do; end do

do i=-1,1; do j=-1,1
 select case(ggd_cell(i+i0,j+j0)%region) 
  case('smth  ')
   uuu(i,j)=uu(i+i0,j+j0)
  case('crit  ')
   ads=ggd_cell(i+i0,j+j0)%ccpt(1)%address
   call visit(ads,temp)
   g_cell=temp%g_cell; p_cell=temp%p_cell
   uuu(i,j)=temp%p_cell%or_state
   if(temp%l_stk%cv_nmb.gt.0.or.temp%r_stk%cv_nmb.gt.0) then
    if(temp%l_stk%cv_nmb.gt.0) then
     ads_n=temp%l_stk
     uuu(i,j)=uuu(i,j)-p_cell%l_state
    else
     ads_n=temp%r_stk
     uuu(i,j)=uuu(i,j)-p_cell%r_state
    end if
    call visit(ads_n,temp_n)
    g_cell_n=temp_n%g_cell; p_cell_n=temp_n%p_cell
    uuu(i,j)=uuu(i,j)+p_cell_n%or_state
   end if
  case('node  ')
   if(i+i0.eq.ndd(nd_nmb)%n_cell%x_idx.and.j+j0.eq.ndd(nd_nmb)% &
      n_cell%y_idx) then
    uuu(i,j)=ndd(nd_nmb)%n_cell%or_state
   else
    uuu(i,j)=ndd(nd_nmb_n)%n_cell%or_state
   end if    
 end select
end do; end do

!conservation=0.d0
!do i=-1,0; do j=-1,0
! conservation=conservation+uuu(i,j)
!end do; end do

conservation=0.d0
do i=-1,1; do j=-1,1
 conservation=conservation+uuu(i,j)
end do; end do 

print*, 'The conserved quantity is', conservation

end subroutine check_conservation


end module node_move_n_split