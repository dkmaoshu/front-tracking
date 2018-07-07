% 	program showc_hr
	load d:\workf9~1\show\showp_hr.dat;
	load d:\workf9~1\show\shown.dat;
	x_width=shown(1);
	y_width=shown(2);
	nxll=shown(3);
	nyll=shown(4);
	nx=shown(5);
	ny=shown(6);
    
    left=10*nxll-5;
    right=10*(nxll+nx-1)+5;
    bottom=10*nyll-5;
    top=10*(nyll+ny-1)+5;
    nnx=10*nx+1;
    nny=10*ny+1;
        
	h=x_width/(nx-1)/10.0;
	A = reshape(showp_hr,nny,nnx);
	for i=1:nnx
	  x(i) = (i+left-1)*h;
	end;
	for i=1:nny
	  y(i) = (i+bottom-1)*h;
 	end;
%	contour(x,y,A,30);
%   axis equal tight
%	axis([-1. 1. -1. 1.]);
%	axis('square');
    mesh(x,y,A);
    view(-35,45);