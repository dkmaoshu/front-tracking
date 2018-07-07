% 	program sect_x_p

    row=6;
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
	  xxx(i) = (i+nxll-1)*h;
	end;
	for i=1:nx
	  yyy(i) = A(row,i);
    end;
    plot(xxx,yyy,'O');
    
    hold on;
    
    clear;
    
    row=6;
    load d:\workf90_1\show1\showu1.dat;
	load d:\workf90_1\show1\shown1.dat;
	x_width1=shown1(1);
	y_width1=shown1(2);
	nxll1=shown1(3);
	nyll1=shown1(4);
	nx1=shown1(5);
	ny1=shown1(6);
	h1=x_width1/(nx1-1);
	B = reshape(showu1,ny1,nx1);
   	for i=1:nx1
	  xxx1(i) = (i+nxll1-1)*h1;
	end;
	for i=1:nx1
	  yyy1(i) = B(row,i);
    end;
    plot(xxx1,yyy1,'black -');    
    
    hold on
    
    xh(1)=-0.124004;
    xh(2)=-0.124004;
    yv(1)=-0.1; 
    yv(2)=0.45

    plot(xh,yv,'R -.');
    
    clear
    
    xh(1)=0.257654;
    xh(2)=0.257654;
    yv(1)=-0.1; 
    yv(2)=0.45;
    
    plot(xh,yv,'R -.');       
    
    axis([-1.2 1.2 -0.1 0.45]);

    hold off