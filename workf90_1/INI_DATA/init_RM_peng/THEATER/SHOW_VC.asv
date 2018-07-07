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
	A = reshape(showv,ny,nx);
	for i=1:nx
	  x(i) = (i+nxll-1)*h;
	end;
	for i=1:ny
	  y(i) = (i+nyll-1)*h;
 	end;

%   B=A((ny-1)/2+70:ny-100,:);
%   yy=y((ny-1)/2+70:ny-100);   
%   contour(x,yy,B,30); 
     
%   contour(x,y,A,30);
   contourf(x,y,A,30,'LineColor','none');
   axis equal tight
 
   axis([0.0,2.0,-2.0,2.0]);

   clear
   
