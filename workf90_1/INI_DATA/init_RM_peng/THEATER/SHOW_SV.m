%%Spike velocity %%%

load d:\workf90_1\output\spike.dat
load d:\workf90_1\output\cur_t.dat

    x1=cur_t;
   
   s1=spike(:,1);
   s2=spike(:,2);
   s3=spike(:,4);

   sv=s1+s3.*(s1-s2)/(-1);
   %sv means Spike velocity

   x1=10*10^(0.5)*x1-25.0;
   sv=100*10^(0.5)*sv;
 
  plot(x1,sv,'b');
  %hold on;

%   axis equal tight
%   axis ([0.0,9.66,-0.5,0.1]);
%   axis ([0.275,0.325,0.45,0.50])
%   axis ([0.60,0.70,0.40,0.50]);
  % axis ([0.0,28.0,-0.06,0.06]);
   %axis ([0.0,800,0.0,18.0]);
   % axis ([0,400,0,18]);
    axis ([0.0,2000.0,0.0,18.0]);
   axis ('square');
%   hold off;
   %hold on;