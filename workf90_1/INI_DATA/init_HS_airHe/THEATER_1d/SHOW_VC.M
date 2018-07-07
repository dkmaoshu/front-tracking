% 	program showv
	load d:\workf90_1\show\showv.dat;
	load d:\workf90_1\show\shown.dat;
	x_width=shown(1);
	y_width=shown(2);
	nxll=shown(3);
	nyll=shown(4);
	nx=shown(5);
	ny=shown(6);
	h=x_width/(nx-1);
    a=0.5*y_width;
    b=0.5*x_width;
	A = reshape(showv,ny,nx);
	for i=1:nx
	  x(i) = (i+nxll-1)*h;
	end;
	for i=1:ny
	  y(i) = (i+nyll-1)*h;
 	end;
	contour(x,y,A,30);
   axis equal tight
%   axis([-2.225 2.225 -0.445 0.445]);
%   axis([0.0,1.0,-0.1,0.9]);
   axis([-b b -a a]);
%   axis([-0.11125,0.0,0.05,0.16125]);
%   axis([-0.2225,-0.11125,-0.11125,0.0]);
%   axis([0.0,0.11125,-0.2225,-0.11125]);
   

   
