%	program showmd

 	load d:\workf9~1\show\shown.dat;
	load d:\workf9~1\show\showm.dat;
	n=shown;
    m=showm;
    nn=n(4)+n(6)-1;
    h=1./n(4);

% 	for i=-nn-1:nn;
% 	  x29(1)=-1.05 ; x29(2)=1.05;
%	  y29(1)=(i+0.5)*h ; y29(2)=(i+0.5)*h;
%	  plot(x29,y29,'r'); hold on;
%	end;
%	for j=-nn-1:nn;
%	  x30(1)=(j+0.5)*h ; x30(2)=(j+0.5)*h;
%	  y30(1)=-1.05 ; y30(2)=1.05;
%	  plot(x30,y30,'r');
%	end;
 
   load d:\workf9~1\show\showmd.dat;
   x3=showmd(:,1);
   y3=showmd(:,2);
   plot(x3,y3,'r');

   axis equal tight;
%   axis([-0.0 1.0 -0.0 1.0]);
%   axis([-0.2225,-0.0,-0.2225,0.0]);
%   axis([0.0,0.445,0.0,0.445]);
%   axis([0.1,0.2225,0.1,0.2225]);  
    axis([-0.445,0.445,-0.445,0.445]);

%	axis('square');
	hold off;
