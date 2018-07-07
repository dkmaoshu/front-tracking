% 	program circulation

    load d:\workf90_1\show\circulation_HS.dat;
    
    x=circulation_HS(:,1);
    y=circulation_HS(:,2);
    plot(x,y,'*');
	
 %   axis([-1.0 1.0 0.9 1.9]);

    hold off