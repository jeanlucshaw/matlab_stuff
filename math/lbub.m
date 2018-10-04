function [lb,ub]	= lbub(data,alpha)


I	= find(isfinite(data)) ;
data	= data(I) ;
N	= numel(data) ;
NI	= floor(N*alpha) ; 
data	= sort(data) ;
ub	= data(end-NI) ;
lb	= data(NI+1) ;

end

