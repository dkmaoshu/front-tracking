% 	program sect_x_p

    row=6;
	load d:\workf90_1\show\showp.dat;
	load d:\workf90_1\show\shown.dat;
	x_width=shown(1);
	y_width=shown(2);
	nxll=shown(3);
	nyll=shown(4);
	nx=shown(5);
	ny=shown(6);
	h=x_width/(nx-1);
	A = reshape(showp,ny,nx);
   	for i=1:nx
	  xxx(i) = (i+nxll-1)*h;
	end;
	for i=1:nx
	  yyy(i) = A(row,i);
    end;
    plot(xxx,yyy,'O');

    hold on
    
%    row=105;
%    load c:\workf90\system\show\showvv.dat;
%	load c:\workf90\system\show\shown.dat;
%	xe_width=shown(1);
%	ye_width=shown(2);
%	nxlle=shown(3);
%	nylle=shown(4);
%	nxe=shown(5);
%	nye=shown(6);
%	hh=x_width/(nxe-1);
%	AA = reshape(showvv,nye,nxe);
%   	for i=1:nxe
%	  xxxe(i) = (i+nxlle-1)*hh;
%	end;
%	for i=1:nye
%	  yyye(i) = AA(row,i);
%    end;
%    plot(xxxe,yyye,'r');
	axis([-1.2 1.2 0.9 1.9]);
%	axis([-0.3 0.3 1.5762 1.6]);
    

    hold off