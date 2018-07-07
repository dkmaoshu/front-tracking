% least square fit
%=====================================
load d:\workf90_1\output\bubble.dat
load d:\workf90_1\output\spike.dat
load d:\workf90_1\output\cur_t.dat

 x1=cur_t;
 y1=bubble(:,1);
 y2=bubble(:,2);
 y3=bubble(:,4);

 z1=spike(:,1);
 z2=spike(:,2);
 z3=spike(:,4);

 yb=y1+y3.*(y1-y2)/(-1);
 zs=z1+z3.*(z1-z2)/(-1);

 %Gr means growth rate
 gr=0.5*(zs-yb);
 x1=10*10^(0.5)*x1-25.0;
 gr=100*10^(0.5)*gr;
 
   m=length (x1);
   x2=x1 (2:m);
%   z1 (1:m)=1;
   x3=x1(1:m-1);
   x4=x2-x3;
   g2=gr(1:m-1);
   g3=gr(2:m);
   a=g2'*x4+g3'*x4;
   q1 (1:m-1)=1.0;
   b=q1*x4;
   c=0.5*a/b;
   
 plot(x1,gr,'b');
 hold on;

 %  Air-SF6  %
 %axis([0.02,30.00,-0.06,0.06]);
 axis([0.0,800.0,0.0,18.0]);
 axis('square');
 hold on;
 