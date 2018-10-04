function y = oreo(x,vp,va)
% Synthax:           y = oreo(x,vp,va)
% 
% Carefree vector concatenation. Prepends 'vp' and appends 'va'
% to vector 'x' and returns in 'y'. If 'x' is a column vector,
% 'y' is a column vector and vice versa.

if iscolumn(x)
	if ~iscolumn(vp); vp = vp'; end
	if ~iscolumn(va); va = va'; end
	y = [vp;x;va] ;
else
	if  iscolumn(vp); vp = vp'; end
	if  iscolumn(va); va = va'; end
	y = [vp x va] ;
end
end
