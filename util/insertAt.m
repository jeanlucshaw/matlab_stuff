function output = insertAt(input,index,value)
%
% output = insertAt(input,index,value)
%
% this function returns the 'input' vector with 'value' inserted
% at position 'index' . Picture it as 'bumping up' the values at
% each of the index positions.
%

N             = numel(input) ;
Ni            = numel(index) ;
Nv            = numel(value) ;

index         = sort(index,'ascend') ;
c             = false(N+Nv*Ni,1);
indexvec      = [] ;
for ii = 1:Ni
	indexvec      = [indexvec index(ii)+(ii-1)*Nv:1:index(ii)+(ii-1)*Nv+Nv-1] ;
end
c(indexvec)   = true;
output        = nan(size(c));
output(~c)    = input;
output(c)     = value ;

end
