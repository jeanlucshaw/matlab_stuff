function qA = boxav(x,y,qx,qy,q)
% Synthax :            qA = boxav(x,y,qx,qy,q)
%
% Box averaging on a rectangular but not necessarily regular grid. Grid
% coordinates are given by matrices 'x' and 'y'. Averaged variable 'q' is
% located at 'qx' and 'qy' coordinates and can be vectors or matrices.
% Averaged variable 'qA', the size of 'x' and 'y' is the output.

% Initialize output
qA = nan(size(x)) ;

% Conduct averaging
for ii = 1:numel(x(:,1))     % longitude loop
    for jj = 1:numel(y(1,:)) % latitude loop
    
    % calculate distance to neighbor points    
    try dxl = abs( x(ii,jj)-x(ii,jj-1) ); catch; dxl = abs( x(ii,jj)-x(ii,jj+1) ); end 
    try dxr = abs( x(ii,jj)-x(ii,jj+1) ); catch; dxr = abs( x(ii,jj)-x(ii,jj-1) ); end 
    try dyu = abs( y(ii,jj)-y(ii+1,jj) ); catch; dyu = abs( y(ii,jj)-y(ii-1,jj) ); end
    try dyd = abs( y(ii,jj)-y(ii-1,jj) ); catch; dyd = abs( y(ii,jj)-y(ii+1,jj) ); end 
    
    % set logical condition to find points of interest
    I = find(qx <= x(ii,jj) + dxr./2 & ... 
             qx >= x(ii,jj) - dxl./2 & ...
             qy <= y(ii,jj) + dyu./2 & ...
             qy >= y(ii,jj) - dyd./2) ;
         
    % average
    if ~isempty(I)
        qA(ii,jj) = mean(q(I),'omitnan') ;
    else
        qA(ii,jj) = nan ;
    end
    
    end   
end

end