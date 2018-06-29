function out = nanmedian(in)

I   = find(isfinite(in)) ;
out = median(in(I)) ;

end
