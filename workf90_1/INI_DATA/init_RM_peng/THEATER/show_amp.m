 
%	program show Amplitude

 load d:\workf90_1\output\amp.dat; 
 load d:\workf90_1\output\cur_t.dat;
	x3=cur_t;
    x3=10*(10^0.5)*x3-25.0;
       
    y3=amp;
    y3=10*y3;
	plot(x3,y3,'r');
	
axis('square')
  
  % axis([0.0,1.0,-1.0,1.0]);
   %  Air-SF6  %
 %axis([00.02,27.90,0.15,1.10]);
 axis([0.0,800.0,0.0,12.0]);
 %axis([-0.5,2000.0,1.5,18.0]);
 %hold on;