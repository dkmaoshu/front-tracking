% 	program showu
	load d:\workf90_1\show\showu.dat;
	load d:\workf90_1\show\shown.dat;
	x_width=shown(1);
	y_width=shown(2);
	nxll=shown(3);
	nyll=shown(4);
	nx=shown(5);
	ny=shown(6);
	h=x_width/(nx-1);
	A = reshape(showu,ny,nx);
	for i=1:nx
	  x(i) = (i+nxll-1)*h;
	end;
	for i=1:ny
	  y(i) = (i+nyll-1)*h;
 	end;
	contourf(x,y,A,30,'LineColor','none');
    axis equal tight
   
   axis([-0.445 0.445 -0.445 0.445]);
%   axis([-2.225,2.225,-0.445,0.445])
%   axis([0.0,1.0,-0.1,0.9]);
%   axis([-3.115,3.115,-0.445,0.445]);
%  axis([0.0,0.11125,-0.2225,-0.11125]);
   
