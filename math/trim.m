function I	= trim(x,p)
% Syntax	I	= trim(x)
%
% Return index values of the top  and bottom 'p'% values in
% 'x'. 'x' may be N dimensional and contain nan values. Nan 
% values are untouched. 'I' is always in linear index format. 

% Vectorize and De-nan
xv	= x(:) ;
xI	= 1:numel(xv) ;
I	= find(isnan(xv)) ;
xv(I)	= [] ;
xI(I)	= [] ;

% Find indexes of values to trim
N	= numel(xv) ;
Ncut	= ceil((p/100)*N) ;
sx	= sort(xv) ;
I	= find(xv < sx(Ncut+1) | xv > sx(N-Ncut)) ; % +1 because indices start at 1
I	= xI(I) ;

end

