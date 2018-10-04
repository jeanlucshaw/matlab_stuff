function y = app(x,v)
% Synthax:             y = app(x,v)
%
% Carefree vector concatenation. 'v' will be appended to 'x'
% along the working dimension and the result is output to 'y'.
% If 'x' is a column vector, 'y' is column vector and vice versa.
if iscolumn(x)
	if ~iscolumn(v); v = v'; end
	y = [x;v] ;
else
	if  iscolumn(v); v = v'; end
	y = [x v] ;
end
end
