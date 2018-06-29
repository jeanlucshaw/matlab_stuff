function time = time_in_seconds(end_time,start_time)
%
%          time = time_in_seconds(end_time,start_time)
%
% Takes as input two time values in datenum format and returns
% the time difference in seconds. end_time must be greater than
%                         start_time .
%

total_raw = end_time - start_time ;
total_str = datestr(total_raw,'yyyy:mm:dd:HH:MM:SS') ;
str_vec   = textscan(total_str,repmat('%f',[1 6]),'delimiter',':') ;
year      = 365*24*3600 ;  % a year in seconds
month     = 30.5*24*3600 ; % an average month in seconds
day       = 24*3600 ;      % a day in seconds
hour      = 3600 ;         % an hour in seconds
minute    = 60 ;           % a minute in seconds
time      = str_vec{1}*year + str_vec{2}*month  + str_vec{3}*day + ...
            str_vec{4}*hour + str_vec{5}*minute + str_vec{6} ;

end
