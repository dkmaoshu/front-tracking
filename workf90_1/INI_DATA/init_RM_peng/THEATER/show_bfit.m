% Best fit curve with spike velocity
 % Best fit curve with spike velocity
%======================================
load d:\workf90_1\output\cur_t.dat

 x1=cur_t;
%  y1=bubble(:,1);
%  y2=bubble(:,2);
%  y3=bubble(:,4);

%  bv = y1+y3.*(y1-y2)/(-1);
 %bv stands for bubble velocity
 

 
  x1=10*10^(0.5)*x1-25.0;
  m=length(x1);
  sv(1:m)=1.0;

  sv = 207./x1.^(0.43);
  %sv = 510./x1.^(0.58);
  % bv=100*10^(0.5)*bv;
 
  plot(real(x1),real(sv),'b--');

 axis([0 2000 0 18]);
 axis ('square');