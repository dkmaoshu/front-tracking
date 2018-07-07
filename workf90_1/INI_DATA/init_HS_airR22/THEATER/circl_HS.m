% 	program circulation

	load d:\workf90_1\show\circulation_HS.dat;
	
    x=circulation_HS(:,1);
    y=circulation_HS(:,2);
    x=x-12.047;
    plot(x,y,'*');

	axis([0.0 260.0 0.0 0.2]);

    hold off