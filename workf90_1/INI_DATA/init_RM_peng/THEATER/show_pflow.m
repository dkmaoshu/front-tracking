% Potential flow model (Hecht et al model)
%===========================================

%load d:\workf90_1\output\bubble.dat
load d:\workf90_1\output\cur_t.dat

 k=-1.675;
 x1=cur_t;%bubble;%
 
 l=length(x1);
 x2=x1(500:l);
 x2=10*10^(0.5)*x2-25.0;
 
%   y1=bubble(:,1);
%   y2=bubble(:,2);
%   y3=bubble(:,4);
% 
%   bv = y1+y3.*(y1-y2)/(-1);
 %bv stands for bubble velocity
 x1=10*10^(0.5)*x1-25.0;
 m=length(x1);
 bv(1:m)=1.0;
 bv=bv';
 %bv=20000*bv./(3*k*x1);
 bv=20000./(3*k*x2);
%  x1(1)
%  pause
%  s=size(x1)
%  pause
%  x1(s(1))
%  pause
  %x1=10*10^(0.5).*x1-35.0;
%   bv=100*10^(0.5)*bv;

  plot(x2,bv,'--');
  
% axis ([0.0,800,-16.0,1.0])
  axis ([0,2000,-18.0,0.0]);
  axis ('square');
 