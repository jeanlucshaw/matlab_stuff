function [yb] = bin_data(x,y,xb)

% Synthaxe : [yb] = bin_data(x,y,xb)

% Cette fonction produit un vecteur 'yb' bined aux bins 'xb'

dx = xb(2) - xb(1) ;

yb    = nan(size(xb));
for i = 1:length(xb)
    j = find(x <= xb(i) + 0.5*dx & x >= xb(i)-0.5*dx) ;
    yb(i) = mean(y(j),'omitnan') ;
end

end
