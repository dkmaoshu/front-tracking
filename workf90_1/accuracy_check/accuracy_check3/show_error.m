% 	program showerror
	load d:\workf90_1\accuracy_check\accuracy_check3\error.dat;
	load d:\workf90_1\accuracy_check\accuracy_check3\shown.dat;
	nxll=shown(1);
	nyll=shown(2);
	nx=shown(3);
	ny=shown(4);
	h=x_width/(nx-1);
	A = reshape(error,ny,nx);
	for i=1:nx
	  x(i) = (i+nxll-1)*h;
	end;
	for i=1:ny
	  y(i) = (i+nyll-1)*h;
 	end;
	mesh(x,y,A,30);
   axis equal tight
	axis([-1.0 1.0 -1.0 1.0]);
