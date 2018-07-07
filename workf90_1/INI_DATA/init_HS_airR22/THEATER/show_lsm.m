% Program growth rate
   load d:\workf90_1\output\initgrs.dat
   load d:\workf90_1\output\initgrb.dat
   load d:\workf90_1\output\initti.dat
   x1=initti
   y1=initgrs(:,1);
   y2=initgrs(:,2);
   y3=initgrs(:,4);
   z1=initgrb(:,1);
   z2=initgrb(:,2);
   z3=initgrb(:,4);
   y4=y3.*(y2-y1)+y1;
   z4=z3.*(z2-z1)+z1;
   g1=0.5*(y4-z4);
   x1=10*10^(0.5)*x1-20.0;
   g1=100*10^(0.5)*g1;
   
   
   m=length(x1);
   x2=x1(2:m);
%   z1(1:m)=1;
   x3=x1(1:m-1);
   x4=x2-x3;
   g2=g1(1:m-1);
   g3=g1(2:m);
   a=g2'*x4+g3'*x4;
   q1(1:m-1)=1.0;
   b=q1*x4;
   c=0.5*a/b;
   
   
   plot(x1,g1,'g');
   hold on;

%   axis equal tight
%   axis([0.0,9.66,-0.5,0.1]);
%   axis([0.275,0.325,0.45,0.50])
%   axis([0.60,0.70,0.40,0.50]);
%   axis([0.0,1.0,-0.1,0.9]);
   axis([0.0,300,-60.0,0.0]);
   axis('square');
%   hold off;
    hold on;