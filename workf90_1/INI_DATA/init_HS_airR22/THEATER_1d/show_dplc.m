%	program showp

 	load d:\workf9~1\show\shown.dat;
	load d:\workf9~1\show\showm.dat;
	n=shown;
        m=showm;
        nn=n(4)+n(6)-1;
        h=1./n(4);

 	for i=-nn-1:nn;
 	  x29(1)=-1.05 ; x29(2)=1.05;
	  y29(1)=(i+0.5)*h ; y29(2)=(i+0.5)*h;
	  plot(x29,y29,'r'); hold on;
	end;
	for j=-nn-1:nn;
	  x30(1)=(j+0.5)*h ; x30(2)=(j+0.5)*h;
	  y30(1)=-1.05 ; y30(2)=1.05;
	  plot(x30,y30,'r');
	end;
 
	load d:\workf9~1\show\showdc1.dat;
	x3=showdc1(:,1);
	y3=showdc1(:,2);
	plot(x3,y3,'b');
	hold on;

   axis equal tight
   axis([-0.3 -0.2 0.45 0.55]);
    
%   pause
   
	if m > 1
	load d:\workf9~1\show\showdc2.dat;
	x4=showdc2(:,1);
	y4=showdc2(:,2);
	plot(x4,y4,'b');
	hold on;
    end

	if m > 2
	load d:\workf9~1\show\showdc3.dat;
	x5=showdc3(:,1);
	y5=showdc3(:,2);
	plot(x5,y5,'b');
        hold on;
	end

	if m > 3
	load d:\workf9~1\show\showdc4.dat;  
	x6=showdc4(:,1);
	y6=showdc4(:,2);
	plot(x6,y6,'b');
        hold on;
	end

	if m > 4
    	load d:\workf9~1\show\showdc5.dat;         
	x9=showdc5(:,1);
	y9=showdc5(:,2);
	plot(x9,y9,'b');
        hold on;
	end


	if m > 5
    	load d:\workf9~1\show\showdc6.dat;         
	x10=showdc6(:,1);
	y10=showdc6(:,2);
	plot(x10,y10,'b');
        hold on;
	end
   
	if m > 6
    	load d:\workf9~1\show\showdc7.dat;         
	x11=showdc7(:,1);
	y11=showdc7(:,2);
	plot(x11,y11,'b');
        hold on;
	end
   
	if m > 7
    	load d:\workf9~1\show\showdc8.dat;         
	x12=showdc8(:,1);
	y12=showdc8(:,2);
	plot(x12,y12,'b');
        hold on;
	end
   
	if m > 8
    	load d:\workf9~1\show\showdc9.dat;         
	x13=showdc9(:,1);
	y13=showdc9(:,2);
	plot(x13,y13,'b');
        hold on;
	end
   
	if m > 9
    	load d:\workf9~1\show\showdc10.dat;         
	x14=showdc10(:,1);
	y14=showdc10(:,2);
	plot(x14,y14,'b');
        hold on;
	end
   
	if m > 10
    	load d:\workf9~1\show\showdc11.dat;         
	x15=showdc11(:,1);
	y15=showdc11(:,2);
	plot(x15,y15,'b');
        hold on;
	end
   
    hold off;
