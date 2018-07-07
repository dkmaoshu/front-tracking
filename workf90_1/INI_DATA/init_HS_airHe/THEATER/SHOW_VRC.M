% 	program show_vrc
   load d:\workf90_1\show\velo_r.dat;
	load d:\workf90_1\show\shown.dat;
	x_width=shown(1);
	y_width=shown(2);
	nxll=shown(3);
	nyll=shown(4);
	nx=shown(5);
	ny=shown(6);
	h=x_width/(nx-1);
	A = reshape(velo_r,ny,nx);
	for i=1:nx
	  x(i) = (i+nxll-1)*h;
	end;
	for i=1:ny
	  y(i) = (i+nyll-1)*h;
 	end;
	contour(x,y,A,30);
   axis equal tight
	axis([-1. 1. -1. 1.]);
