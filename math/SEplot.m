function [bins,bs,m,M]=SEplot(dist)
% Synthax :            [bins,bs,m,M]=SEplot(dist)
%     
% Produces a histogram plot of distribution 'dist' with a color code
% identifying which portions of the distribution correspond to +-1 standard
% error, +-2 standard error, +-3 standard error, and beyond. Limits are
% included to the categories and standard error is calculated with the
% built in matlab std function

sigma = std(dist,'omitnan') ;
mdist = mean(dist,'omitnan') ;

M     = max(dist) ; 
m     = min(dist) ;

bs    = sigma./20 ;
bins  = sort([mdist:-bs:m mdist+bs:bs:M]);

figure ; 
i = find(dist <= mdist + 1*sigma & dist >= mdist - 1.*sigma) ;
histogram(dist(i),bins,'facecolor','b') ;

hold on

i = find(dist <= mdist + 2*sigma & dist >  mdist + 1.*sigma | ... 
         dist <  mdist - 1*sigma & dist >= mdist - 2.*sigma) ;
histogram(dist(i),bins,'facecolor','g') ;

i = find(dist <= mdist + 3*sigma & dist >  mdist + 2.*sigma | ... 
         dist <  mdist - 2*sigma & dist >= mdist - 3.*sigma) ;
histogram(dist(i),bins,'facecolor','r') ;

i = find(dist <= mdist - 3*sigma | dist >= mdist + 3.*sigma) ;
histogram(dist(i),bins,'facecolor','c') ;

end
