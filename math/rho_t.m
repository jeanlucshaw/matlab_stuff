function [rho] = rho_t(x,y,degrad)
%                      [rho] = rho_t(x,y,degrad)
%
% Explicitly calculates the correlation coefficient rho_t detailed in
% Fisher and Lee (1983), between the two circular variables 'x' and 'y'.
% Variable 'degrad' is a string and should be "deg" if 'x' and 'y' are
% in degrees or 'rad' if they are in radians.
%
% Reference: Fisher, N.I., LEE, A.J.,(1983), A correlation coefficient
% 		for circular data, Biometrika, 70, 2, pp. 327-32


numerator	 = 0 ;
denA		 = 0 ;
denB		 = 0 ;

switch degrad
	case 'deg'
		x = x*pi/180 ;
		y = -y*pi/180 ;
	case 'rad'
		% do nothing
end

for ii = 1:numel(x)-1
	numerator = numerator + sum( sin(x(ii)-x(ii+1:end)).*sin(y(ii)-y(ii+1:end)) ) ;
	denA	  = denA      + sum( sin(x(ii)-x(ii+1:end)).^2 );
	denB	  = denB      + sum( sin(y(ii)-y(ii+1:end)).^2 ) ;
end

rho 	= numerator./( sqrt(denA).*sqrt(denB) ) ;

end
