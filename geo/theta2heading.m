
function output  = theta2heading(input)
% Synthax :           output = theta2heading(input)
%
% Takes as input angle vector in cartesian coordinates running counter-
% clockwise from 0 (east) to 360 degrees and transforms it into a heading
% vector running clockwise from 0 to the north to 360 degrees. 
%
 
input            = input - 90;
input            = 360 - input ;
input(input < 0) = input(input < 0) + 360 ;
output            = mod(input,360) ;

end
