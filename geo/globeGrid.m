function [LNG,LAT] = globeGrid(lims,sz)
% Synthax :         [LNG,LAT] = globeGrid(lims,sz)
%
% Takes as input longitude and latitude limits 'lims', in format :
%
%             lims = [lon_min lon_max lat_min lat_max]
% 
% and a grid size 'sz' in kilometers, and generates a grid uniform in
% kilometer space, with the center of the x dimension aligned with the
% center of the area specified by 'lims'. This is mostly designed for use
% with small grids (~50-100km^2) where this representation looks almost
% regular in coordinate space and kind of "looks nice"... at least to me it
% does.
%
% 

% EARTH RADIUS
R       = 6371 ;

% convert to km
lng = [lims(1); lims(1); lims(2); lims(2)] ;
lat = [lims(3); lims(4); lims(3); lims(4)] ;
avl = mean(lng) ;
lng = lng - avl ;                              % move to prime meridian

% SINUSOIDAL PROJECTION
y     = (            lat*R*pi / 180 ) ;
x     = ( cosd(lat).*lng*R*pi / 180 ) ;

avy   = mean(y) ;
avx   = mean(x) ;

% MAKE GRID VECTORS AROUND CENTER IN KM SPACE
xg    = [flipud( [avx:-sz:min(x)]' ) ; [avx+sz:sz:max(x)]' ] ;
yg    = [flipud( [avy:-sz:min(y)]' ) ; [avy+sz:sz:max(y)]' ] ;

% GENERATE UNIFORM GRID IN KM SPACE
[X,Y] = meshgrid(xg,yg) ;

% CONVERT BACK TO GEOGRAPHICAL COORDINATES
LAT   = 180*Y./(R*pi) ;
LNG   = 180*X./(R*pi*cosd(LAT)) + avl ;   % moved back to average meridian

end