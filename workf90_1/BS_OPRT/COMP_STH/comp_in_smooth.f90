module computation_in_smooth
! This module computes the numerical solution in smooth regions, 
! which includes solution in ordinary grid cells and the left and 
! right states in critical cells.

use solu_comput
! 'solu_comp.f90'

use numerical_fluxes
! 'flux.f90'

use node_cells
! 'nd_cls.f90

use boundary_conditions
! 'boundary_conditions.f90'

implicit none

type crit_info
 type(adss_info) :: address
 type(geo_info) :: g_cell
 type(phy_info) :: p_cell
 type(state), dimension(0:2) :: l_state, r_state, or_state
 character*6 :: side
end type crit_info

type(crit_info) :: crit

type(state), dimension(:,:), allocatable :: flx, fly, flxm, flym
type(state), dimension(:,:), allocatable :: uum, uu1

public  compute_smooth
private uu, acvv, ggd_cell, smth_dpr, compute_in_smooth, &
        compute_on_curves, compute_flux_coef, pick_smth_flux, &
        coeff, data_fill, flx, fly, flxm, flym, uum, uu1, &
		comput_ndfl1, crit, artificial_pressure, boundary_condition1, &
		boundary_condition2


contains


subroutine compute_smooth

implicit none

integer :: i, j
   
allocate(uu(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(acvv(curves_number)) 
allocate(ggd_cell(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(ndd(ndn))

allocate(flx(nxll-2:nxll+nx,nyll-2:nyll+ny))
allocate(fly(nxll-2:nxll+nx,nyll-2:nyll+ny))
allocate(flxm(nxll-2:nxll+nx,nyll-2:nyll+ny))
allocate(flym(nxll-2:nxll+nx,nyll-2:nyll+ny))
allocate(uum(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))
allocate(uu1(nxll-3:nxll+nx+2,nyll-3:nyll+ny+2))

!$omp parallel do
do i=nxll-2,nxll+nx
 do j=nyll-2,nyll+ny
  flx(i,j)=error_data; fly(i,j)=error_data
  flx(i,j)%gamma=error_data; fly(i,j)%gamma=error_data
  flxm(i,j)=error_data; flym(i,j)=error_data
  flxm(i,j)%gamma=error_data; flym(i,j)%gamma=error_data
 end do
end do
!$omp end parallel do

!$omp parallel do
do i=nxll-3,nxll+nx+2
 do j=nyll-3,nyll+ny+2 
  uum(i,j)=error_data; uu1(i,j)=error_data
  uum(i,j)%gamma=error_data; uu1(i,j)%gamma=error_data
 end do
end do
!$omp end parallel do

call give_sth(0)
call give_acv
call give_gcl
call give_nd

! call scanning(uu)

! call check_list_ac(1,'down  ')

call boundary_condition1

! ggd_cl=ggd_cell(-1,-3)
! call check_list_ac(1,'down  ')

call compute_in_smooth

! call check_list_ac(1,'down  ')

call compute_on_curves(1)

if(boundary_type_2.ne.'none        ') call nonreflecting_boundary_1(uum)

! call scanning(uu)

! call check_list_ac(1,'down  ')

uu1=uu; uu=uum

!$omp parallel do
do i=nxll-3,nxll+nx+2
 do j=nyll-3,nyll+ny+2
  uum(i,j)=error_data
 end do
end do
!$omp end parallel do

flxm=flx; flym=fly

! call check_list_ac(1,'down  ')

call boundary_condition2

call compute_in_smooth

! call check_list_ac(1,'down  ')

call compute_on_curves(2)

if(boundary_type_2.ne.'none        ') call nonreflecting_boundary_1(uum)

! call check_list_ac(1,'down  ')

! call scanning(uum)

!$omp parallel do
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  uu(i,j)=0.5d0*(uu1(i,j)+uum(i,j))
 end do
end do
!$omp end parallel do

!$omp parallel do
do i=nxll-2,nxll+nx
 do j=nyll-2,nyll+ny
  flx(i,j)=0.5d0*(flx(i,j)+flxm(i,j))
  fly(i,j)=0.5d0*(fly(i,j)+flym(i,j))
 end do
end do
!$omp end parallel do

!call hidden_state_check_a

!call artificial_pressure

! call check_list_ac(1,'down  ')

call boundary_condition2

if(boundary_type_2.ne.'none        ') call nonreflecting_boundary_2(uum,uu1)

! call check_list_ac(1,'down  ')

call compute_flux_coef

! call scanning(uu1)

call comput_ndfl1

! call check_list_ac(1,'down  ')

call get_sth(1)
call get_acv
call get_nd
call get_gcl


! ccc=crit

deallocate(uu); deallocate(acvv)
deallocate(ggd_cell); deallocate(ndd)

deallocate(flx); deallocate(fly); deallocate(flxm)
deallocate(flym); deallocate(uum); deallocate(uu1)

end subroutine compute_smooth


subroutine compute_in_smooth
! This subroutine updates the numerical solution in ordinary
! grid cells.

implicit none

integer	:: i, j
type(auxi_crit_cell), pointer :: temp
!type(crit_info_a) :: crit
type(state), dimension(-3:3) :: ux, uy
type(state), dimension(-2:3) :: data_sector
real(8) :: x_posi, y_posi

!$omp parallel do private(temp,ux,uy,data_sector) 
do i=nxll,nxll+nx-1
 do j=nyll,nyll+ny-1
  if(ggd_cell(i,j)%region.eq.'smth  ') then
!   call clean_up_g_cell(crit%g_cell); call clean_up_p_cell(crit%p_cell)
!   crit%address=adss_info(0,0); crit%side=' ' 
   nullify(temp)
   call smth_dpr(i,j,temp,ux,'xx ','      ','ordinary',.true.,3,'forward ')
   call smth_dpr(i,j,temp,uy,'yy ','      ','ordinary',.true.,3,'forward ')
! Produce the data-sectors needed in evaluation of numerical fluxes.
   x_posi=(dfloat(i)+0.5d0)*h; y_posi=dfloat(j)*h
   data_sector=ux(-2:3); flx(i,j)=flux(data_sector,'xx',x_posi,y_posi)
   x_posi=dfloat(i)*h; y_posi=(dfloat(j)+0.5d0)*h
   data_sector=uy(-2:3); fly(i,j)=flux(data_sector,'yy',x_posi,y_posi)
   call smth_dpr(i,j,temp,ux,'xx ','      ','ordinary',.true.,3,'backward')
   call smth_dpr(i,j,temp,uy,'yy ','      ','ordinary',.true.,3,'backward')
! Produce the data-sectors needed in evaluation of numerical fluxes.
   x_posi=(dfloat(i)-0.5d0)*h; y_posi=dfloat(j)*h
   data_sector=ux(-3:2); flx(i-1,j)=flux(data_sector,'xx',x_posi,y_posi)
   x_posi=dfloat(i)*h; y_posi=(dfloat(j)-0.5d0)*h
   data_sector=uy(-3:2); fly(i,j-1)=flux(data_sector,'yy',x_posi,y_posi)
! Evaluate numerical fluxes.
   uum(i,j)=uu(i,j)-r*(flx(i,j)-flx(i-1,j)+fly(i,j)-fly(i,j-1))
!   uum(i,j)%gamma=uu(i,j)%gamma
  end if
 end do
end do
!$omp end parallel do

end subroutine compute_in_smooth


subroutine compute_on_curves(level)
! This subroutine updates the left and right states on 
! discontinuity curves.

implicit none
integer, intent(in) :: level

integer	:: i, j
!type(crit_info_a) :: crit
type(state), dimension(-3:3) :: ux, uy
type(state), dimension(-2:3) :: data_sector
integer :: k, kk 
type(auxi_crit_cell), pointer :: temp, boundary_end
integer :: head_mark, end_mark
logical :: head_switch
type(state) :: uc, uc1
type(state), dimension(4) :: fl
type(geo_info) :: g_cell
real(8) :: x_posi, y_posi

! type(auxi_crit_cell), pointer :: txxx

! txxx=>acvv(1)%eend%previous

do k=1,curves_number
 if(acvv(k)%status.eq.'awake') then
  if(associated(acvv(k)%begin%next)) then
   temp=>acvv(k)%begin%next
   head_mark=temp%address%idx; head_switch=.true.
   boundary_end=>acvv(k)%eend_boundary%previous
   if(acvv(k)%eend%previous%address.eq.boundary_end%address) then
    end_mark=error_index
   else
    end_mark=acvv(k)%eend%previous%next%address%idx
   end if
    
   do while(associated(temp).and.head_switch)
    
    g_cell=temp%g_cell
   	  
    uc=temp%p_cell%l_state
    call smth_dpr(i,j,temp,ux,'xx ','left  ','critical',.true.,3,'forward ')
    call smth_dpr(i,j,temp,uy,'yy ','left  ','critical',.true.,3,'forward ')
! Compute the data-sectors needed in evaluation of numerical fluxes.
    x_posi=(dfloat(g_cell%x_idx)+0.5d0)*h
    y_posi=dfloat(g_cell%y_idx)*h
    data_sector=ux(-2:3); fl(2)=flux(data_sector,'xx',x_posi,y_posi)
	x_posi=dfloat(g_cell%x_idx)*h
	y_posi=(dfloat(g_cell%y_idx)+0.5d0)*h
    data_sector=uy(-2:3); fl(3)=flux(data_sector,'yy',x_posi,y_posi)
    call smth_dpr(i,j,temp,ux,'xx ','left  ','critical',.true.,3,'backward')
    call smth_dpr(i,j,temp,uy,'yy ','left  ','critical',.true.,3,'backward')
! Compute the data-sectors needed in evaluation of numerical fluxes.
    x_posi=(dfloat(g_cell%x_idx)-0.5d0)*h
	y_posi=dfloat(g_cell%y_idx)*h
    data_sector=ux(-3:2); fl(4)=flux(data_sector,'xx',x_posi,y_posi)
	x_posi=dfloat(g_cell%x_idx)*h
	y_posi=(dfloat(g_cell%y_idx)-0.5d0)*h
    data_sector=uy(-3:2); fl(1)=flux(data_sector,'yy',x_posi,y_posi)
    uc1=uc-r*(fl(2)-fl(4)+fl(3)-fl(1))
! Evaluate numerical fluxes.
    temp%c_cell%l_state(level)=uc1
    if(level.eq.1) then
     temp%f_cell%l_flux=fl
    else
     temp%c_cell%l_state(2)=0.5d0*(temp%c_cell%l_state(2)+ &
     temp%c_cell%l_state(0))
	 do kk=1,4
      temp%f_cell%l_flux(kk)=temp%f_cell%l_flux(kk)+fl(kk)
	  temp%f_cell%l_flux(kk)=0.5d0*temp%f_cell%l_flux(kk)
     end do
    end if 
! Compute left states on discontinuity curves.
    
!    crit%side='right'
    uc=temp%p_cell%r_state
    call smth_dpr(i,j,temp,ux,'xx ','right ','critical',.true.,3,'forward ')
    call smth_dpr(i,j,temp,uy,'yy ','right ','critical',.true.,3,'forward ')
! Compute the data-sectors needed in evaluation of numerical fluxes.
    x_posi=(dfloat(g_cell%x_idx)+0.5d0)*h
	y_posi=dfloat(g_cell%y_idx)*h
    data_sector=ux(-2:3); fl(2)=flux(data_sector,'xx',x_posi,y_posi)
	x_posi=dfloat(g_cell%x_idx)*h
	y_posi=(dfloat(g_cell%y_idx)+0.5d0)*h
    data_sector=uy(-2:3); fl(3)=flux(data_sector,'yy',x_posi,y_posi)
    call smth_dpr(i,j,temp,ux,'xx ','right ','critical',.true.,3,'backward')
    call smth_dpr(i,j,temp,uy,'yy ','right ','critical',.true.,3,'backward')
! Compute the data-sectors needed in evaluation of numerical fluxes.
    x_posi=(dfloat(g_cell%x_idx)-0.5d0)*h
	y_posi=dfloat(g_cell%y_idx)*h
    data_sector=ux(-3:2); fl(4)=flux(data_sector,'xx',x_posi,y_posi)
	x_posi=dfloat(g_cell%x_idx)*h
	y_posi=(dfloat(g_cell%y_idx)-0.5d0)*h
    data_sector=uy(-3:2); fl(1)=flux(data_sector,'yy',x_posi,y_posi)
! Evaluate numerical fluxes.
    uc1=uc-r*(fl(2)-fl(4)+fl(3)-fl(1))
    temp%c_cell%r_state(level)=uc1
    if(level.eq.1) then
     temp%f_cell%r_flux=fl
    else
     temp%c_cell%r_state(2)=0.5d0*(temp%c_cell%r_state(2)+ &
     temp%c_cell%r_state(0))
	 do kk=1,4
      temp%f_cell%r_flux(kk)=temp%f_cell%r_flux(kk)+fl(kk)
	  temp%f_cell%r_flux(kk)=0.5d0*temp%f_cell%r_flux(kk)
     end do
    end if 
! Compute right states on discontinuity curves.    
    
    temp=>temp%next
    if(associated(temp)) then
     head_switch=((temp%address%idx.ne.head_mark).and.(temp%address%idx.ne.end_mark))
    else
     exit
    end if
   end do
  end if
 end if
end do
      
do k=1,curves_number
 if(acvv(k)%status.eq.'awake') then
  if(associated(acvv(k)%begin%next)) then
   temp=>acvv(k)%begin%next
   head_mark=temp%address%idx; head_switch=.true.
   boundary_end=>acvv(k)%eend_boundary%previous
   if(acvv(k)%eend%previous%address.eq.boundary_end%address) then
    end_mark=error_index
   else
    end_mark=acvv(k)%eend%previous%next%address%idx
   end if
   do while(associated(temp).and.head_switch)
	temp%p_cell%l_state=temp%c_cell%l_state(level)
    temp%p_cell%r_state=temp%c_cell%r_state(level)
    temp=>temp%next
    if(associated(temp)) then
     head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
    else
     exit
    end if
   end do
  end if
 end if

!  call check_list_ac(5,'up    ')

end do

end subroutine compute_on_curves


subroutine artificial_pressure
! Smooth parts of the solution may collapse near the tracked discontinuities. Artificial 
! pressure is used to prevent the collapse.

implicit none

integer	:: i, j, k
type(geo_info) :: g_cell, g_cell_p
type(phy_info) :: p_cell
type(comp_info) :: c_cell
type(state), dimension(-3:3) :: uxy
type(auxi_crit_cell), pointer :: temp, boundary_end
integer :: head_mark, end_mark
logical :: head_switch
character*8 :: partner

do k=1,curves_number

 if(acvv(k)%status.eq.'awake') then
  if(associated(acvv(k)%begin%next)) then
   temp=>acvv(k)%begin%next
   head_mark=temp%address%idx; head_switch=.true.
   boundary_end=>acvv(k)%eend_boundary%previous
   if(acvv(k)%eend%previous%address.eq.boundary_end%address) then
    end_mark=error_index
   else
    end_mark=acvv(k)%eend%previous%next%address%idx
   end if
   do while(associated(temp).and.head_switch)
	g_cell=temp%g_cell; c_cell=temp%c_cell
	c_cell%left_pressure=error_data
	c_cell%right_pressure=error_data
	temp%c_cell=c_cell
	temp%c_cell=c_cell
    temp=>temp%next
    if(associated(temp)) then
     head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
    else
     exit
    end if
   end do
  end if
 end if

 if(acvv(k)%status.eq.'awake') then
  if(associated(acvv(k)%begin%next)) then
   temp=>acvv(k)%begin%next
   head_mark=temp%address%idx; head_switch=.true.
   boundary_end=>acvv(k)%eend_boundary%previous
   if(acvv(k)%eend%previous%address.eq.boundary_end%address) then
    end_mark=error_index
   else
    end_mark=acvv(k)%eend%previous%next%address%idx
   end if
   do while(associated(temp).and.head_switch)

	g_cell=temp%g_cell; p_cell=temp%p_cell; c_cell=temp%c_cell
	partner=temp%partner; call clean_up_g_cell(g_cell_p)
	if(partner.eq.'previous') then
	 g_cell_p=temp%previous%g_cell
    else
	 g_cell_p=temp%next%g_cell
    end if

! Artificially pressed in $x$- and $y$-direction.
	call compressed('left  ')
	call compressed('right ')
     
!	call compressed('left  ','xx ')
!    call compressed('left  ','yy ')	    
!    call compressed('right ','xx ')
!    call compressed('right ','yy ')
	  
	temp%c_cell=c_cell
    temp=>temp%next
    if(associated(temp)) then
     head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
    else
     exit
    end if

   end do
   
  end if
 end if

 if(acvv(k)%status.eq.'awake') then
  if(associated(acvv(k)%begin%next)) then
   temp=>acvv(k)%begin%next
   head_mark=temp%address%idx; head_switch=.true.
   boundary_end=>acvv(k)%eend_boundary%previous
   if(acvv(k)%eend%previous%address.eq.boundary_end%address) then
    end_mark=error_index
   else
    end_mark=acvv(k)%eend%previous%next%address%idx
   end if
   do while(associated(temp).and.head_switch)

	g_cell=temp%g_cell; p_cell=temp%p_cell; c_cell=temp%c_cell
	if(c_cell%left_pressure.gt.0.9d0*error_data) then
	 p_cell%l_state=p_cell%l_state+c_cell%left_pressure
    end if
	if(c_cell%right_pressure.gt.0.9d0*error_data) then
	 p_cell%r_state=p_cell%r_state+c_cell%right_pressure
    end if

	temp%p_cell=p_cell
    temp=>temp%next
    if(associated(temp)) then
     head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
    else
     exit
    end if

   end do
   
  end if
 end if

end do


contains


!subroutine compressed(side,direction)
subroutine compressed(side)

implicit none
character*6, intent(in) :: side
!character*3, intent(in) :: direction

!logical :: if_needed
character*8 :: for_backward_x, for_backward_y
integer :: single
character*6 :: other_side
type(state) :: uc, uc1, aa_x, aa_y, aa
type(state), dimension(4) :: data_sector_x, data_sector_y
!real(8) :: dd
type(geo_info) :: g_cell_4_check
real(8) :: alpha, beta, rr
integer :: edge_number
real(8), dimension(2) :: pt1, pt2

! Preparation.
select case(side)
 case('left  '); uc=p_cell%l_state
 case('right '); uc=p_cell%r_state
 case default; call error_message
end select
other_side=side_shift(side)
!if_needed=.false.

! First, determine whether an artificial pressure in this direction is needed.
!if(g_cell%g_type.eq.direction) then
! if_needed=.true.
!else
! if(g_cell%g_type.eq.'xy ') if_needed=.true.
! if(g_cell%g_type.eq.'xy ') then
!  call find_single(g_cell,single)
!  if(g_cell%point(single).eq.side) if_needed=.true.
! end if 
!end if

! Second, determine the 'for_backward' fashion.
 !if(if_needed) then
! Determine the 'g_type' and the 'points' of the critical cell for check first.
g_cell_4_check=g_cell
if(partner.ne.'single  ') then
 g_cell_4_check%g_type='xy '
 select case(partner)
  case('previous'); edge_number=2
  case('next    '); edge_number=1
 end select
 select case(g_cell%g_type)
  case('xx ')
   if(g_cell_p%x_idx.lt.g_cell%x_idx) then
    select case(g_cell%edge(edge_number))
     case(1); g_cell_4_check%point(4)=side_shift(g_cell%point(4))
     case(3); g_cell_4_check%point(1)=side_shift(g_cell%point(1))
	 case default; call error_message 
    end select
   else
    select case(g_cell%edge(edge_number))
     case(1); g_cell_4_check%point(3)=side_shift(g_cell%point(3))
     case(3); g_cell_4_check%point(2)=side_shift(g_cell%point(2))
	 case default; call error_message
    end select
   end if
  case('yy ')
   if(g_cell_p%y_idx.lt.g_cell%y_idx) then
    select case(g_cell%edge(edge_number))
     case(2); g_cell_4_check%point(1)=side_shift(g_cell%point(1))
     case(4); g_cell_4_check%point(2)=side_shift(g_cell%point(2))
	 case default; call error_message
    end select
   else
    select case(g_cell%edge(edge_number))
     case(2); g_cell_4_check%point(4)=side_shift(g_cell%point(4))
     case(4); g_cell_4_check%point(3)=side_shift(g_cell%point(3))
	 case default; call error_message
    end select
   end if
  case('xy '); call error_message   
 end select
end if

! Then determine for_backward_x and for_backward_y.
for_backward_x='        '; for_backward_y='        '
select case(g_cell_4_check%g_type)
! case('xx ','yy ')
 case('xx ')
  if(g_cell_4_check%point(1).eq.side) then
   for_backward_x='backward'
  else
   if(g_cell_4_check%point(1).eq.other_side) then
    for_backward_x='forward '
   else
    call error_message
   end if
  end if
 case('yy ')
  if(g_cell_4_check%point(1).eq.side) then
   for_backward_y='backward'
  else
   if(g_cell_4_check%point(1).eq.other_side) then
    for_backward_y='forward '
   else
    call error_message
   end if
  end if
 case('xy ')
!  if(g_cell_4_check%point(1).eq.side.and.g_cell_4_check%point(4).eq.side) then
!   for_backward_x='backward'
!  else
!   if(g_cell_4_check%point(2).eq.side.and.g_cell_4_check%point(3).eq.side) then
!    for_backward_x='forward '
!   else
!    call error_message
!   end if
!  end if
!  if(g_cell_4_check%point(1).eq.side.and.g_cell_4_check%point(2).eq.side) then
!   for_backward_y='backward'
!  else
!   if(g_cell_4_check%point(3).eq.side.and.g_cell_4_check%point(4).eq.side) then     
!    for_backward_y='forward '
!   else
!    call error_message
!   end if
!  end if	   
  call find_single(g_cell_4_check,single)
  select case(single)
   case(1)
    if(g_cell_4_check%point(single).eq.side) then
     for_backward_x='backward'; for_backward_y='backward'
    else
     for_backward_x='forward '; for_backward_y='forward '
    end if
   case(2)
    if(g_cell%point(single).eq.side) then
     for_backward_x='forward '; for_backward_y='backward'
    else
     for_backward_x='backward'; for_backward_y='forward '
    end if
   case(3)
    if(g_cell%point(single).eq.side) then
     for_backward_x='forward '; for_backward_y='forward '
    else
     for_backward_x='backward'; for_backward_y='backward'
    end if
   case(4)
    if(g_cell%point(single).eq.side) then
     for_backward_x='backward'; for_backward_y='forward '
    else
     for_backward_x='forward '; for_backward_y='backward'
    end if
   end select
end select
!end if

! Finally, artificially compress the solution.
!if(if_needed) then 	    
call smth_dpr(i,j,temp,uxy,'xx ',side,'critical',.true.,3,'full    ')
select case(for_backward_x)
 case('forward '); data_sector_x=uxy(0:3)
 case('backward'); do i=0,3; data_sector_x(i+1)=uxy(-i); end do
 case default; do i=1,4; data_sector_x(i)=error_data; end do
end select
aa_x=3.0d0*data_sector_x(2)-3.0d0*data_sector_x(3)+data_sector_x(4)

! call design_pressure(data_sector,uc1)

call smth_dpr(i,j,temp,uxy,'yy ',side,'critical',.true.,3,'full    ')
select case(for_backward_y)
 case('forward '); data_sector_y=uxy(0:3)
 case('backward'); do i=0,3; data_sector_y(i+1)=uxy(-i); end do
 case default; do i=1,4; data_sector_y(i)=error_data; end do
end select
aa_y=3.0d0*data_sector_y(2)-3.0d0*data_sector_y(3)+data_sector_y(4)
  
call pick_point(g_cell,1,pt1); call pick_point(g_cell,2,pt2)
alpha=pt2(2)-pt1(2); beta=pt2(1)-pt1(1)
rr=dsqrt(alpha*alpha+beta*beta)
if(rr.gt.0.1d-6) then
 alpha=dabs(alpha)/rr; beta=dabs(beta)/rr
else
 alpha=0.5d0*dsqrt(2.0d0); beta=0.5d0*dsqrt(2.0d0)
end if

! Compute the artificial pressure.
if(aa_x.gt.0.9d0*error_data.and.aa_y.gt.0.9d0*error_data) then
 aa=alpha*alpha*aa_x+beta*beta*aa_y
else
 if(aa_x.gt.error_data) then
  aa=aa_x
 else
  if(aa_y.gt.error_data) then
   aa=aa_y
  else
   call error_message
  end if
 end if
end if
if(data_sector_x(1).gt.0.9d0*error_data) then
 uc1=(aa-data_sector_x(1))/8.0d0
else
 if(data_sector_y(1).gt.0.9d0*error_data) then
  uc1=(aa-data_sector_y(1))/8.0d0
 else
  call error_message
 end if
end if      
! dd=dabs(aa)
! dd=dd/h
! select case(g_cell%g_type)
!  case('xx ','yy '); dd=dmin1(dd*dd,1.0d0)
!  case('xy ');       dd=dmin1(dd*dd,1.0d0)
!  case default; call error_message
! end select

! select case(for_backward)
!  case('forward '); uc1=(-r)*dd*(uxy(0)-uxy(1))
!  case('backward'); uc1=(-r)*dd*(uxy(0)-uxy(-1))
!  case default; call error_message
! end select

select case(side)
 case('left  ')
! if(c_cell%left_pressure.lt.0.9d0*error_data) then
  c_cell%left_pressure=uc1
! else
!  c_cell%left_pressure=(c_cell%left_pressure+uc1)/2.0d0
! end if
 case('right ')
!  if(c_cell%right_pressure.lt.0.9d0*error_data) then
   c_cell%right_pressure=uc1
!  else
!   c_cell%right_pressure=(c_cell%right_pressure+uc1)/2.0d0
!  end if
end select

!end if 

end subroutine compressed


end subroutine artificial_pressure


subroutine comput_ndfl1
! This subroutine computes numerical fluxes on cell-edges of node
! cells that lies in smooth region.

implicit none
type(node_info) :: n_cell
integer :: i, j, edge_nb, ii
integer, dimension(2) :: ind 


do i=1, ndn
 if(ndd(i)%status.eq.'awake ') then
  n_cell=ndd(i)%n_cell
  do ii=1,4; n_cell%flux(ii)=error_data; end do

  do j=1, 4
! Find the indexes of the adjacent cell.
   select case(j) 
    case(1)
     ind(1)=n_cell%x_idx
     ind(2)=n_cell%y_idx-1
     edge_nb=3
    case(2)  
     ind(1)=n_cell%x_idx+1
     ind(2)=n_cell%y_idx
     edge_nb=4
    case(3)
     ind(1)=n_cell%x_idx 
     ind(2)=n_cell%y_idx+1
     edge_nb=1
    case(4)
     ind(1)=n_cell%x_idx-1 
     ind(2)=n_cell%y_idx
     edge_nb=2
   end select

   if(ggd_cell(ind(1),ind(2))%region.eq.'smth  ') then
! If the adjacent cell is located in the smooth region, pick 
! the numerical flux across the cell-edge.
    select case(edge_nb)
     case(1); n_cell%flux(j)=fly(ind(1),ind(2)-1)
     case(2); n_cell%flux(j)=flx(ind(1),ind(2))
     case(3); n_cell%flux(j)=fly(ind(1),ind(2))
     case(4); n_cell%flux(j)=flx(ind(1)-1,ind(2))
    end select
   end if

  end do

  ndd(i)%n_cell=n_cell
 end if
end do

end subroutine comput_ndfl1


subroutine compute_flux_coef
! In computation in critical cells, numerical fluxs on partial 
! cell-edges are required; thus, polynomials of second order are
! required to describe the numerical fluxes on cell-edges. This 
! subroutine computes coefficients of these polynomials.
implicit none

integer	:: k, kk, kx
type(auxi_crit_cell), pointer :: temp, boundary_end
type(flux_info) :: f_cell
type(state), dimension(-3:3) :: flxy
!real(8), dimension(-3:3) :: flxy_t
type(state), dimension(0:2) :: cff, cfb
type(state), dimension(-1:1) :: fxy
integer :: head_mark, end_mark
logical :: head_switch

 integer :: iii, jjj

do k=1,curves_number
 if(acvv(k)%status.eq.'awake') then
  if(associated(acvv(k)%begin%next)) then
   temp=>acvv(k)%begin%next
   head_mark=temp%address%idx; head_switch=.true.
   boundary_end=>acvv(k)%eend_boundary%previous
   if(acvv(k)%eend%previous%address.eq.boundary_end%address) then
    end_mark=error_index
   else
    end_mark=acvv(k)%eend%previous%next%address%idx
   end if

   do while(associated(temp).and.head_switch)
    crit%g_cell=temp%g_cell; crit%p_cell=temp%p_cell
    crit%address=temp%address; f_cell=temp%f_cell

    crit%side='left'
    do kk=1,4
     call pick_smth_flux(f_cell,flxy,kk)
     call data_fill_flux(flxy,3)
!	 do kx=-3,3; flxy_t(kx)=flxy(kx)%value; end do
!	 call data_fill(flxy_t,3)
!	 do kx=-3,3; flxy(kx)%value=flxy_t(kx); end do
	 do kx=-1,1; fxy(kx)=flxy(kx); end do
     call coeff(fxy,cff,cfb)
     do kx=0,2
	  f_cell%l_cff(kk,kx)=cff(kx)
	  f_cell%l_cfb(kk,kx)=cfb(kx)
     end do
    end do

     iii=crit%g_cell%x_idx; jjj=crit%g_cell%y_idx

    crit%side='right'
    do kk=1,4
     call pick_smth_flux(f_cell,flxy,kk)
!	 flxy_t=flxy%value
	 call data_fill_flux(flxy,3)
!	 flxy%value=flxy_t
	 do kx=-1,1; fxy(kx)=flxy(kx); end do
     call coeff(fxy,cff,cfb)
	 do kx=0,2
      f_cell%r_cff(kk,kx)=cff(kx)
	  f_cell%r_cfb(kk,kx)=cfb(kx)
     end do
    end do

    temp%f_cell=f_cell

    temp=>temp%next
    if(associated(temp)) then
     head_switch=(temp%address%idx.ne.head_mark.and.(temp%address%idx.ne.end_mark))
    else
     exit
    end if
   end do
  end if
 end if
end do

end subroutine compute_flux_coef


subroutine pick_smth_flux(f_cell,flxy,edge)
! Pick up numerical fluxes on cell-edges that belong to smooth
! region.

implicit none

type(flux_info), intent(in) :: f_cell
! The auxiliary critical cell to work on.
type(state), dimension(-3:3), intent(out) :: flxy
! The data sector to store the picked-up fluxes.
integer, intent(in) :: edge
! The cell-edge to work on.

integer :: i0, j0, k, edge1, i
character*8 :: p_ext, n_ext

do i=-3,3; flxy(i)=error_data; end do
if(crit%side.eq.'left') then
 flxy(0)=f_cell%l_flux(edge)
else
 flxy(0)=f_cell%r_flux(edge)
end if
! Start with the cell-edge under concern.

i0=crit%g_cell%x_idx; j0=crit%g_cell%y_idx
p_ext='denied  '; n_ext='denied  '
if(ggd_cell(i0,j0)%ccpt(edge)%address.eq.crit%address.and. &
crit%g_cell%point(edge).eq.crit%side) then
 if(edge.eq.1.or.edge.eq.2) then
  n_ext='allowed '
 else
  p_ext='allowed '
 end if
end if
edge1=cycle_4(edge+1)
if(ggd_cell(i0,j0)%ccpt(edge1)%address.eq.crit%address.and. &
crit%g_cell%point(edge1).eq.crit%side) then
 if(edge.eq.1.or.edge.eq.2) then
  p_ext='allowed '
 else
  n_ext='allowed '
 end if
end if
! Variables 'p_ext' and 'n_ext' indicate whether the pick-up are
! allowed to be carried on in the positive and negative directions.

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! At present we do not pick the nearby numerical fluxes; therefore, the flux coefficients 
! are calculated in first-order fashion/
 p_ext='denied '; n_ext='denied '

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

if(p_ext.eq.'allowed') then
 select case(edge)
  case(1)
   do k=1,3
    if(fly(i0+k,j0-1).lt.0.9d0*error_data) exit
	flxy(k)=fly(i0+k,j0-1)
   end do
  case(2)
   do k=1,3
    if(flx(i0,j0+k).lt.0.9d0*error_data) exit 
	flxy(k)=flx(i0,j0+k)
   end do
  case(3)
   do k=1,3
    if(fly(i0+k,j0).lt.0.9d0*error_data) exit
	flxy(k)=fly(i0+k,j0)
   end do
  case(4)
   do k=1,3
    if(flx(i0-1,j0+k).lt.0.9d0*error_data) exit
	flxy(k)=flx(i0-1,j0+k)
   end do
 end select
end if
! Pick up numerical fluxes in the left smooth region for the 
! four cell edges.

if(n_ext.eq.'allowed') then
 select case(edge)
  case(1)
   do k=-1,-3, -1
    if(fly(i0+k,j0-1).lt.0.9d0*error_data) exit
	flxy(k)=fly(i0+k,j0-1)
   end do
  case(2)
   do k=-1,-3, -1
    if(flx(i0,j0+k).lt.0.9d0*error_data) exit
	flxy(k)=flx(i0,j0+k)
   end do
  case(3)
   do k=-1,-3, -1
    if(fly(i0+k,j0).lt.0.9d0*error_data) exit 
	flxy(k)=fly(i0+k,j0)
   end do
  case(4)
   do k=-1,-3, -1
    if(flx(i0-1,j0+k).lt.0.9d0*error_data) exit
	flxy(k)=flx(i0-1,j0+k)
   end do
 end select
end if
! Pick up numerical fluxes in the right smooth region for the 
! four cell-edges.

end subroutine pick_smth_flux


subroutine coeff(fxy,cff,cfb)
! The computation of the coefficients of the polynomials for
! numerical fluxes.
implicit none
type(state), intent(in) :: fxy(-1:1)
type(state), intent(out) :: cff(0:2), cfb(0:2)
cff(0)=(5.d0*fxy(0)-fxy(1))/8.d0; 
cff(1)=fxy(0); cff(2)=(fxy(1)-fxy(0))/2.d0
cfb(0)=(3.d0*fxy(0)+fxy(-1))/8.d0 
cfb(1)=fxy(0); cfb(2)=(fxy(0)-fxy(-1))/2.d0
end subroutine coeff


end module computation_in_smooth