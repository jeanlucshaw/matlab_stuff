function output = nanmean(input,varargin)

%
% Matlab is super cheap for not including this 
%     in the student license. Seriously. 
%

I = find(isfinite(input)) ;
output = mean(input(I),varargin{:}) ;
end
