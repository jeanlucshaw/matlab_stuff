function [lond,latd] = dm2deci(londeg,lonmin,lonsec,EW,latdeg,latmin,latsec,NS)

%
% This function converts coordinates from deg/min/sec to decimal
%
% [lond,latd] = dm2deci(londeg,lonmin,lonsec,EW,latdeg,latmin,latsec,NS)
%
% Enter numbers for longitude degrees, minutes and seconds, and a string for cardinal direction ('E' or 'W')
%                                              followed by
%       numbers for latitude  degrees, minutes and seconds, and a string for cardinal direction ('N' or 'S')
%
% Written by Jean-Luc Shaw
%

if(strmatch(NS,'N')) ; lat_sign =  1. ; end
if(strmatch(NS,'S')) ; lat_sign = -1. ; end
if(strmatch(EW,'E')) ; lon_sign =  1. ; end
if(strmatch(EW,'W')) ; lon_sign = -1. ; end

lond = lon_sign*(londeg+lonmin./60+lonsec./3600) ;
latd = lat_sign*(latdeg+latmin./60+latsec./3600) ;

end
