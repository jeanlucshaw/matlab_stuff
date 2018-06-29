function [t_out,idx,num] = fillGaps(t_in,step,tol)
%
% Synthax :      [t_out,idx,num] = fillGaps(t_in,step,tol)
%
% Utility designed to find gaps in reasonably regular time series and
% insert a patch of regularly grid times where the hole was. 't_out' will
% therefore be a longer version of the vector 't_in' where a patch given
% by:
%
%               [t_in(ii) + step :step: t_in(ii+1) - step]
% 
% or it's transpose, with a gap greater than:
%
%                             step*(1+tol)
%
% having been found between t_in(ii) and t_in(ii+1), has been wedged
% between those two values in the time series. The output values 'idx' and
% 'num' respectively indicate the index where you need to insert a patch
% into a vector the orginial size of t_in, and how big a patch, so that it
% then still matches t_out. For example by doing:
%
% for ii = 1:numel(idx); h = insertAt(h,idx(ii),nan([num(ii) 1])) ; end
%

idx    = [] ;
num    = [] ;

% FIND HOLES
I      = find( diff(t_in) > step*(2+tol) ) ;

for ii = 1:numel(I) 
    
    % GENERATE AND INSERT PATCH
    pat        = [t_in(I(ii))+step:step:t_in(I(ii)+1)-step] ;
    if ~isempty(pat)
        if iscolumn(t_in) ; pat = pat' ; end
        t_in       = insertAt(t_in,I(ii)+1,pat) ;

        % KEEP IN MEMORY WHERE IT WENT AND HOW BIG IT WAS
        idx        = [idx; I(ii)+1] ;
        num        = [num; numel(pat)] ;

        % UPDATE LOCATION OF HOLES TO PATCH
        I          = I + num(end) ;
    end
    
end

% YOU ARE DONE PATCHING, OUTPUT THE SERIES
t_out = t_in ;

end