function theta = vvangle(v1,v2)
% Synthax :                theta = vvangle(v1,v2)
%
% Get the absolute angle between two line vectors, or between the vectors
% described by the lines v1 and v2.
%

%% make sure the sizes are compatible
if ~( isequal(size(v1),size(v2)) || isvector(v1) && isvector(v2) && numel(v1) == numel(v2) )
    disp('    Error : v1 and v2 have different dimensions!') ;
    return ;
end

%% get the vector norms
nv1 = sum(abs(v1).^2,2).^(1/2) ;
nv2 = sum(abs(v2).^2,2).^(1/2) ;

%% calculate the angle
theta = acosd( ( dot(v1,v2,2) ) ./ (nv1.*nv2) ) ;

end