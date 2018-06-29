function [dist,I] = SEtrunc(dist,thres)
% Synthax :            [dist,I] = SEtrunc(dist,thres)
% 
% Replaces with nan the values in dist which are outside the range
% 
%                      +-'thres'*std('dist')
%
% where std is matlab built in function. Limits of the range are excluded.
% also returns the index values 'I' of the values nan-ed.

sigma = std(dist,'omitnan') ;
mdist = mean(dist,'omitnan') ;

I       = find(dist < mdist - thres.*sigma | dist > mdist + thres.*sigma) ;
dist(I) = nan ;

end


