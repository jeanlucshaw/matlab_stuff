function [u,v,spd,head] = llt2spd(lon,lat,time)

% Synthax : [u,v,spd,head] = llt2spd(lon,lat,time)
% 
% Returns the speed 'spd' and it's eastward/northward components 'u' and
% 'v' in (m/s), and the heading 'head' in degrees with 0 as east increasing
% towards the north (90 degrees). 
%
% Values are calculated from a latitude/longitude time series with 
% coordinates in decimal degrees and time in days.
% 
% Values are estimated between longitude/latitude coordinates and then
% interpolated linearly to the input 'time' vector. Begining and end
% values are extrapolated using the "pchip" method.

%% make everything column vectors
if ~iscolumn(lon)  ; lon  = lon'  ; end
if ~iscolumn(lat)  ; lat  = lat'  ; end
if ~iscolumn(time) ; time = time' ; end

%% make differential vectors
dlon = diff(lon) ;
dlat = diff(lat) ;
dt   = diff(time) ;
dl   = m_lldist(lon,lat)*1000 ; % (m)

%% calculate norm of speed
spd = dl./(dt*24*60*60) ;

%% calculate heading
head = range_pm180_2_360(atan2d(dlat,dlon)) ;

%% decompose into eastward/northward
u = spd.*cosd(head) ;
v = spd.*sind(head) ;

%% make the time vector
t = time(1:end-1)+0.5*dt(1);

%% interpolate to the input time grid
u    = interp1(t,u,time,'linear','extrap') ;
v    = interp1(t,v,time,'linear','extrap') ;
spd  = interp1(t,spd,time) ;
spd(1)   = sqrt(u(1).^2+v(1).^2);
spd(end) = sqrt(u(end).^2+v(end).^2);
head     = interp1(t,head,time) ;


end
