%%Bubble best fit line

% load d:\workf90_1\output\bubble.dat
% load d:\workf90_1\output\spike.dat
load d:\workf90_1\output\cur_t.dat

 x1=cur_t;
%  y1=bubble(:,1);
%  y2=bubble(:,2);
%  y3=bubble(:,4);

%  z1=spike(:,1);
%  z2=spike(:,2);
%  z3=spike(:,4);

%  yb=y1+y3.*(y1-y2)/(-1);
%  zs=z1+z3.*(z1-z2)/(-1);

 %Gr means growth rate
 %gr=0.5*(zs-yb);
 x1=10*10^(0.5)*x1-25.0;
%  yb=100*10^(0.5)*yb;
 m=length(x1);
 yb(1:m)=1.0;
 yb=4370./x1.^(0.8135);
 plot(real(x1),real(yb),'r');

 %   axis equal tight
 
  % Air_He
%   axis ([0.0,9.66,-0.5,0.1]);
%   axis ([0.275,0.325,0.45,0.50])
%   axis ([0.60,0.70,0.40,0.50]);
  % axis ([0.0,28.0,-0.06,0.06]);
  
  %%%% Air_SF6  %%%%%%
   %axis ([0.0,450,-16.0,0.00]);
   axis ([0,2000,-16.0, 0]);
   axis ('square');
%   hold off;
   %hold on;