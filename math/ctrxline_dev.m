clear all;

m = 1.1;
b = 0 ;

x1 = [0:0.1:5] ;
y1 = m*x1+b ;

x2 = 3+sind([0:5:300]+50) ;
y2 = 3+cosd([0:5:300]+50) ;

%plot(x1,y1,'r') ;
hold on;
plot(x2,y2,'g.')

[Xout,Yout]	= ctrxline(m,b,x2,y2) ;

plot(Xout,Yout,'+') ;

return;

% INITIALISE OUTPUT
Xout = [] ;
Yout = [] ;

% LOOP OVER CONTOUR VERTICES

	% Make the parametric equation for this line
	vec	= @(x)m*x + b ;

	% Draw for the next contours x values
	yvec	= vec(x2) ;
plot(x2,yvec) ;

	% Find intersections
	I	= false(size(x2)) ;
	for jj = 2:numel(x2)
		if y2(jj) > yvec(jj) & y2(jj-1) <= yvec(jj-1) | y2(jj) < yvec(jj) & y2(jj-1) >= yvec(jj-1)
			I(jj) = true ;
		end
	end

	% Boundary
	if y2(1) > yvec(1) & y2(end) <= yvec(end) | y2(1) < yvec(1) & y2(end) >= yvec(end)
		I(1) = true ;
	end

	I	= find(I) ;

	% SOLVE FOR THE INTERCEPTS
	for jj = 1:numel(I) 
		if I(jj) == 1
			P	= polyfit([x2(I(jj)) x2(end)],[y2(I(jj)) y2(end)],1) ;
		else
			P	= polyfit([x2(I(jj)) x2(I(jj)-1)],[y2(I(jj)) y2(I(jj)-1)],1) ;
		end

		% Solve the linear system of these two lines
		A 	= [ m -1;  P(1) -1];
		B	= [-b	; -P(2)	];	
		X	= A\B ;
		Xout	= app(Xout,X(1)) ;
		Yout	= app(Yout,X(2)) ;
	end


plot(Xout,Yout,'+')
