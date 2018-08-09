function [u,v] = ctrspd(x1,y1,x2,y2,t) 
% Synthax :         [u,v] = ctrspd(x1,y1,x2,y2,t)
%
% Calculates the speed at which a contour travels perpendicularly to itself.
% The instantaneous speed of a vertex is calculated as the distance between
% itself at time t_1, and the intersection of the contour at time t_2 with a
% line normal to contour 1 and crossing the vertex.
%
% 	Inputs are coordinates at t_1 : x1,y1
%		   coordinates at t_2 : x2,y2
%      		   times 1 & 2        : t 
%
% Where t is size two vector. Contours may be periodic or non periodic. The
% beginning and ends of the contour will be handled differently in these
% cases.

% PARAMETERS 
NA	= numel(x1) ;
NB	= numel(x2) ;
per1	= false;
per2	= false;
warning('off','MATLAB:singularMatrix') ;
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale') ;
warning('off','MATLAB:nearlySingularMatrix') ;

% MAKE CONTOURS PERIODIC according to prior periodicity
if     x1(1) == x1(end) & y1(1) == y1(end) & x2(1) == x2(end) & y2(1) == y2(end) 
	x1 = [x1(end-1) x1(1:end-1) x1(1)] ;
	y1 = [y1(end-1) y1(1:end-1) y1(1)] ;
	x2 = [x2(end-1) x2(1:end-1) x2(1)] ;
	y2 = [y2(end-1) y2(1:end-1) y2(1)] ;
	per1	= true;
	per2	= true;
elseif x1(1) == x1(end) & y1(1) == y1(end) & x2(1) ~= x2(end) & y2(1) ~= y2(end) 
	x1 = [x1(end-1) x1(1:end-1) x1(1)] ;
	y1 = [y1(end-1) y1(1:end-1) y1(1)] ;
	x2 = [x2(end) x2 x2(1)] ;
	y2 = [y2(end) y2 y2(1)] ;
	per1	= true;
elseif x1(1) ~= x1(end) & y1(1) ~= y1(end) & x2(1) ~= x2(end) & y2(1) ~= y2(end) 
	x1 = [x1(end) x1 x1(1)] ;
	y1 = [y1(end) y1 y1(1)] ;
	x2 = [x2(end-1) x2(1:end-1) x2(1)] ;
	y2 = [y2(end-1) y2(1:end-1) y2(1)] ;
	per2	= true;
else
	x1 = [x1(end) x1 x1(1)] ;
	y1 = [y1(end) y1 y1(1)] ;
	x2 = [x2(end) x2 x2(1)] ;
	y2 = [y2(end) y2 y2(1)] ;
end

% OUTPUT VARIABLES
u	= nan(size(x1)) ;
v	= nan(size(x1)) ;

% LOOP OVER CONTOUR VERTICES
for ii = 2:numel(x1) - 1 
	gam = range_pm180_2_360( atan2d(y1(ii+1)-y1(ii-1),x1(ii+1)-x1(ii-1))) ;
	if     ~per1 & ii == 2
		gam = range_pm180_2_360( atan2d(y1(ii+1)-y1(ii),x1(ii+1)-x1(ii))) ;
	elseif ~per1 & ii == numel(x1) - 1
		gam = range_pm180_2_360( atan2d(y1(ii)-y1(ii-1),x1(ii)-x1(ii-1))) ;
	end
	bet = gam - 90 ;
	if bet < 0 ; bet+360 ; end

	xp = x1(ii) + cosd(bet) ;
	yp = y1(ii) + sind(bet) ;

	% Make the parametric equation for this line
	m1	= polyfit([x1(ii) xp],[y1(ii) yp],1) ;
	vec	= @(x)m1(1)*x + m1(2) ;

	% Draw for the next contours x values
	yvec	= vec(x2) ;

	% Find intersections
	I	= false(NB,1) ;
	for jj = 2:numel(x2) - 1
		if y2(jj) > yvec(jj) & y2(jj-1) < yvec(jj-1) | y2(jj) < yvec(jj) & y2(jj-1) > yvec(jj-1)
			I(jj) = true ;
		end
	end
	I	= find(I) ;

	% Narrow it down to the nearest intersection
	d	= nan(size(I)) ;
	for jj = 1:numel(I)
		d(jj) = sqrt( (x1(ii) - x2(I(jj))).^2 + (y1(ii) - y2(I(jj))).^2 ) ;
	end
	[~,Imin] = min(d) ;
	Imin 	 = I(Imin) ;

	% Fit this line segment
	m2	= polyfit([x2(Imin) x2(Imin-1)],[y2(Imin) y2(Imin-1)],1) ;

	% Solve the linear system of these two lines
	A 	= [ m1(1) -1; m2(1) -1];
	B	= [-m1(2); -m2(2)];	
	X	= A\B ;

	% Calculate angle and magnitude of velocity
	dist	= sqrt( (x1(ii) - X(1)).^2 + (y1(ii) - X(2)).^2 ) ;
	the	= range_pm180_2_360( atan2d(X(2)-y1(ii),X(1)-x1(ii))) ;

	% Decompose into components
	spd	= dist/(t(2) - t(1)) ;
	u(ii)	= spd*cosd(the) ;
	v(ii)	= spd*sind(the) ;

end

% DROP PERIODICITY
u([1 end]) = [];
v([1 end]) = [];
x1([1 end]) = [];
y1([1 end]) = [];

% MIRROR INPUT PERIODICITY
if per1 ; u = app(u,u(1)) ; v = app(v,v(1)) ; end  

function y = app(x,v)

if iscolumn(x)
	y = [x;v] ;
else
	y = [x v] ;
end
end

end
