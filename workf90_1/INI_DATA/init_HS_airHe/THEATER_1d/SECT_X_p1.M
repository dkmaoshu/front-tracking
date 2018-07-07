% 	program sect_x_p1

    row=6;
    load d:\workf90_1\show1\showp1.dat;
	load d:\workf90_1\show1\shown1.dat;
	x_width1=shown1(1);
	y_width1=shown1(2);
	nxll1=shown1(3);
	nyll1=shown1(4);
	nx1=shown1(5);
	ny1=shown1(6);
	h1=x_width1/(nx1-1);
	B = reshape(showp1,ny1,nx1);
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
    yv(1)=0.9; 
    yv(2)=1.8;
    
    plot(xh,yv,'R -.');
    
    xh(1)=0.257654;
    xh(2)=0.257654;
    yv(1)=0.9; 
    yv(2)=1.8;
    
    plot(xh,yv,'R -.');   