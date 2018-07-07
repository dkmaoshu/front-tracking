%	program showdp

   load d:\workf90_1\show\shown.dat;
   load d:\workf90_1\show\showm.dat;
   n=shown;
   m=showm;
   nx=n(3)+n(5)-1;
   ny=n(4)+n(6)-1;
   h=2./(n(5)-1);
   
   
  for i=n(4):ny
 	  x29(1)=-2.0 ; x29(2)=2.0;
	  y29(1)=(i+0.5)*h ; y29(2)=(i+0.5)*h;
      plot(x29,y29,'r');
      hold on;      
  end;
  for j=n(3):nx
     x30(1)=(j+0.5)*h ; x30(2)=(j+0.5)*h;%
	 y30(1)=-7.0 ; y30(2)=7.0; hold on;
	 plot(x30,y30,'r');
  end;
 
	load d:\workf90_1\show\showdc1.dat;
	x3=showdc1(:,1);
	y3=showdc1(:,2);
	plot(x3,y3,'b');
	hold on;

   axis equal tight
 
   axis([0.0 2.0 -2.0 2.0]);
%   axis([0.0 2.00 -7.0 7.0]);

   hold off;