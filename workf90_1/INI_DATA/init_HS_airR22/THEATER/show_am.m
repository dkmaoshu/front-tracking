% Program amplitude
   load d:\workf90_1\output\initam.dat
   load d:\workf90_1\output\initti.dat
   x1=initti;
   x1=10*10^(0.5)*x1-20.0;
   y1=initam;
   y1=10*y1;
   plot(x1,y1,'g');
   hold on;

%   axis equal tight
   axis([0.0,300.0,-9.0,2.0]);
%   axis([0.275,0.325,0.45,0.50])
%   axis([0.60,0.70,0.40,0.50]);
%   axis([0.0,1.0,-0.1,0.9]);
   axis('square');
   hold on;