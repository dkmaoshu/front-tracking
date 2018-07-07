% 	program showrho
	load d:\workf90_1\show\showvort.dat;
	load d:\workf90_1\show\shown.dat;
	x_width=shown(1);
	y_width=shown(2);
	nxll=shown(3);
	nyll=shown(4);
	nx=shown(5);
	ny=shown(6);
	h=x_width/(nx-1);
	A = reshape(showvort,ny,nx);
    for k=1:21;
        vv(k)=(k-11)*8.2+0.1;
    end;        
	for i=1:nx
	  x(i) = (i+nxll-1)*h;
	end;
	for i=1:ny
	  y(i) = (i+nyll-1)*h;
 	end;
   contour(x,y,A,30);
   axis equal tight
%   axis([-2.225 2.225 -0.445 0.445]);
   axis([-0.445 0.445 -0.445 0.445]);
%   axis([0.5,0.75,0.25,0.5]);
%   axis([0.0,1.0,-0.1,0.9]);
