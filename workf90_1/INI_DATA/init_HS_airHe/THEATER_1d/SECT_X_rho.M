% 	program sect_x_rho

    row=6;
	load d:\workf90_1\show\showrho.dat;
	load d:\workf90_1\show\shown.dat;
	x_width=shown(1);
	y_width=shown(2);
	nxll=shown(3);
	nyll=shown(4);
	nx=shown(5);
	ny=shown(6);
	h=x_width/(nx-1);
	A = reshape(showrho,ny,nx);
   	for i=1:nx
	  xxx(i) = (i+nxll-1)*h;
	end;
	for i=1:nx
	  yyy(i) = A(row,i);
    end;
    plot(xxx,yyy,'O');

    hold on
    
    clear
    
    row=6;
	load d:\workf90_1\show1\showrho1.dat;
	load d:\workf90_1\show1\shown1.dat;
	x_width=shown1(1);
	y_width=shown1(2);
	nxll=shown1(3);
	nyll=shown1(4);
	nx=shown1(5);
	ny=shown1(6);
	h=x_width/(nx-1);
	B = reshape(showrho1,ny,nx);
   	for i=1:nx
	  xxx(i) = (i+nxll-1)*h;
	end;
	for i=1:nx
	  yyy(i) = B(row,i);
    end;
    plot(xxx,yyy,'black -');    
    
	axis([-1.2 1.2 0.0 1.8]);

    hold off