function A = ll2area(lng,lat,unit)
% Synthax :                A = ll2area(lng,lat,unit)
%
% Takes as input the geographical coordinate vectors for longitude 'lng'
% and latitude 'lat' in column format and returns the area inside the
% non-intersecting polygon they form as ouput 'A' in units 'unit'. 
% recognized unit options are 'm' for square meters and 'km' for square
% kilometers.
%
% BASED ON:
%
% https://www.periscopedata.com/blog/polygon-area-from-latitude-and-longitude-using-sql
%
% and
%
% http://geomalgorithms.com/a01-_area.html

% PARAMETERS
switch unit
	case 'km'
		R   =  6371 ;    % earth radius (km)
	case 'm'
		R   =  6371000 ; % earth radius (m)
end

% CONVERT TO XY
y     = (            lat*R*pi / 180 ) ;
x     = ( cosd(lat).*lng*R*pi / 180 ) ;
avx   = mean(x) ;
avy   = mean(y) ;
normx = x - avx ;
normy = y - avy ;

% GET ANGLE FOR EVERY VERTEX (0-360) and convert 
ang            = atan2d(normy,normx); % + (-1*normx./abs(normx))*90 + 180 ;
ang( ang < 0 ) = ang( ang < 0 ) + 360 ;
ang            = ang - 90;
ang            = 360 - ang ;
ang(ang < 0)   = ang(ang < 0) + 360 ;
ang            = mod(ang,360) ;

% SORT TO AVOID SELF INTERSECTION
[~,IA] = sort(ang,'descend') ;
x      = x(IA) ;
y      = y(IA) ;

% APPEND FOR WRAPPING
x1 = x(1); xn = x(end) ;
y1 = y(1); yn = y(end) ;
x  = [xn; x; x1] ;
y  = [yn; y; y1] ;

% COMPUTE THE SUM
A      = 0 ;
for ii = 2:numel(x)-1
    A = A + 0.5*( x(ii).*( y(ii+1) - y(ii-1) ) ) ;
end

end
