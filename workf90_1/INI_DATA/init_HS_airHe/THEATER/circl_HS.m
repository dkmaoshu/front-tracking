% 	program circulation

	load d:\workf90_1\show\circulation_HS.dat;
	
    x=circulation_HS(:,1);
    y=circulation_HS(:,2);
    x=x-12.055;
    plot(x,y,'*');

%	axis([0.0 100.0 -0.19 0.0]);

    hold off