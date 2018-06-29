function output = smooth2D(input,hs,vs)
% Synthax :          output = smooth2D(input,hs,vs)
%
% 2D moving average. The window size is 'hs' in direction j and 'vs' in
% direction i. Both 'hs' and 'vs' MUST be odd integers. Each element of
% 'output' is calculated from the nanmean of a box of size 'hs' by 'vs' 
% centered on the element of corresponding indices in 'input'.

output = repmat(nan,size(input)) ;
hpad   = (hs-1)/2 ;
vpad   = (vs-1)/2 ;
if vpad < 1 ; vpad = 0 ; end;
if hpad < 1 ; hpad = 0 ; end;
if iseven(vpad) ; vpad = vpad + 1 ; end
if iseven(hpad) ; hpad = hpad + 1 ; end

[II,JJ] = size(input); 

for ii = 1:II
    for jj = 1:JJ
        if     ii <= vpad    & jj > hpad & jj <= JJ-hpad  % TOP 
            box           = input(1:ii+vpad,jj-hpad:jj+hpad) ;
        elseif ii >  II-vpad & jj > hpad & jj <= JJ-hpad  % BOTTOM
            box           = input(ii-vpad:II,jj-hpad:jj+hpad) ;
        elseif jj <= hpad    & ii > vpad & ii <= II-vpad  % LEFT   
            box           = input(ii-vpad:ii+vpad,1:jj+hpad) ;
        elseif jj >  JJ-hpad & ii > vpad & ii <= II-vpad  % RIGHT   
            box           = input(ii-vpad:ii+vpad,jj-hpad:JJ) ;
        elseif ii <= vpad    & jj <= hpad                 % TOP LEFT
            box           = input(1:ii+vpad,1:jj+hpad) ;
        elseif ii >  II-vpad & jj <= hpad                 % BOTTOM LEFT
            box           = input(ii-vpad:II,1:jj+hpad) ;
        elseif ii <= vpad    & jj >  JJ-hpad              % TOP RIGHT
            box           = input(1:ii+vpad,jj-hpad:JJ) ;
        elseif ii >  II-vpad & jj >  JJ-hpad              % BOTTOM RIGHT
            box           = input(ii-vpad:II,jj-hpad:JJ) ;
        else                                % INSIDE
            box           = input(ii-vpad:ii+vpad,jj-hpad:jj+hpad) ;
        end
        output(ii,jj) = mean(box(:),'omitnan') ;
    end
end
end
