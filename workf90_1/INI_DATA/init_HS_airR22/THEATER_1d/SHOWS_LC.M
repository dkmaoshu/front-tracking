% 	program shows_local
	load d:\workf9~1\show\showu.dat;
	load d:\workf9~1\show\shown.dat;
	x_width=shown(1);
	y_width=shown(2);
	nxll=shown(3);
	nyll=shown(4);
	nx=shown(5);
	ny=shown(6);
	h=x_width/(nx-1);
	A = reshape(showu,ny,nx);
    i0=28;
    j0=17;
	for i=1,21;
	  x(i) = i-10;
	end;
	for j=1, 21;
	  y(j) = j-10;
 	end;
    B=A(i0-10:i0+10,j0-10:j0+10);
%    B=fliplr(B);
%    B=flipud(B);
%	A = reshape(showu,nx,ny);
%	axes('position',rect);
%	colormap([1 1 1]);
	mesh(B);
	view(45,45);
