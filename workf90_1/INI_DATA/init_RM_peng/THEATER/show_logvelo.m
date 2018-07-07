%%log velocity %%%

load d:\workf90_1\output\spike.dat
load d:\workf90_1\output\bubble.dat
load d:\workf90_1\output\cur_t.dat

    x1=cur_t;
   
   s1=spike(:,1);
   s2=spike(:,2);
   s3=spike(:,4);

   sv=s1+s3.*(s1-s2)/(-1);
   %sv means Spike velocity

   x1=10*10^(0.5)*x1-25.0;
   sv=100*10^(0.5)*sv;
  x1=log(x1);
  sv=log(sv);
  plot(x1,sv,'b');
  hold on;

%bubble velocity
 y1=bubble(:,1);
 y2=bubble(:,2);
 y3=bubble(:,4);
 yb=y1+y3.*(y1-y2)/(-1); 
 
  yb=100*10^(0.5)*yb;
  yb=log(yb);
  plot(x1,yb,'r--');

  % axis ([0,400,0,18]);
    axis ([5.5,7.8,0.0,3.0]);
   axis ('square');
%   hold off;
   %hold on;