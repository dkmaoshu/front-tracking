module Haas_Sturtevant_boundary

use solu_comput
! 'solu_comp.f90'

implicit none

public  Haas_Sturtevant_1, Haas_Sturtevant_2, Haas_Sturtevant_3
private smooth_region
 

contains


subroutine Haas_Sturtevant_1

implicit none

call smooth_region

end subroutine Haas_Sturtevant_1


subroutine Haas_Sturtevant_2

implicit none

call smooth_region

end subroutine Haas_Sturtevant_2


subroutine Haas_Sturtevant_3

implicit none

end subroutine Haas_Sturtevant_3


subroutine smooth_region

implicit none

integer :: i, j

nxll_boundary=nxll-3; nx_boundary=nx+3
nyll_boundary=nyll-3; ny_boundary=ny+3

do i=nxll-1,nxll-3,-1
 do j=nyll,nyll+ny-1
  if(ggd_cell(i+nx-1,j)%region.eq.'smth  ') then
   uu(i,j)=uu(nxll,j)
   call clean(ggd_cell(i,j))
  end if
 end do
end do

do i=nxll+nx,nxll+nx+2
 do j=nyll,nyll+ny-1
  if(ggd_cell(i-nx+1,j)%region.eq.'smth  ') then
   uu(i,j)=uu(nxll+nx-1,j)
   call clean(ggd_cell(i,j))
  end if
 end do
end do

do i=nxll, nxll+nx-1
 do j=nyll-3, nyll-1
  uu(i,j)=uu(i,j+ny-1)
  call clean(ggd_cell(i,j))
 end do
 do j=nyll+ny,nyll+ny+2
  uu(i,j)=uu(i,j-ny+1)
  call clean(ggd_cell(i,j))
 end do
end do

end subroutine smooth_region


end module Haas_Sturtevant_boundary