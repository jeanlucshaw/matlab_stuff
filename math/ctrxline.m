function [Xout,Yout] = ctrxline(m,h,b,Xin,Yin) 
% Syntax :                   [Xout,Yout] = ctrxline(m,h,b,Xin,Yin) 
%
% Find all intersections between a line defined by 
%
%		y = m(x - h) + b
%
% and an arbitrary contour define the by the vertices
% in 'Xin' and 'Yin'.
%

% INITIALISE OUTPUT
Xout = [] ;
Yout = [] ;

% Make the parametric equation for this line
vec	= @(x)m*(x-h) + b ;

% Draw for the next contours x values
yvec	= vec(Xin) ;

% Find intersections
I	= false(size(Xin)) ;
for jj = 2:numel(Xin)
	if Yin(jj) > yvec(jj) & Yin(jj-1) <= yvec(jj-1) | Yin(jj) < yvec(jj) & Yin(jj-1) >= yvec(jj-1)
		I(jj) = true ;
	end
end

% Boundary
if Yin(1) > yvec(1) & Yin(end) <= yvec(end) | Yin(1) < yvec(1) & Yin(end) >= yvec(end)
	I(1) = true ;
end

I	= find(I) ;

% Solve for the intercepts
for jj = 1:numel(I) 
	% Line equation of the intercepts
	if I(jj) == 1
		P	= polyfit([Xin(I(jj)) Xin(end)],[Yin(I(jj)) Yin(end)],1) ;
	else
		P	= polyfit([Xin(I(jj)) Xin(I(jj)-1)],[Yin(I(jj)) Yin(I(jj)-1)],1) ;
	end

	% Solve the linear system of these two lines
	A 	= [ m -1	;  P(1) -1];
	B	= [-b+m*h	; -P(2)	];	
	X	= A\B ;
	Xout	= app(Xout,X(1)) ;
	Yout	= app(Yout,X(2)) ;

end
end
