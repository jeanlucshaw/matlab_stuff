function S = uvz2S(u,v,z)
% Synthax:            S = uvz2S(u,v,z)
% 
% Wrapper for the gradient function which takes as input
% horizontal velocities 'u' and 'v' as well as their 
% corresponding depths 'z'. Note that 'u','v', and 'z'
% are all assumed to be of same size.
%
% Returns the horizontal shear:
%
%      S = sqrt( (du/dz)^2 + (dv/dz)^2 )
%
D = ndims(u) ;
M = nan([1 D]) ;
M = size(u) ;
I = find( M > 1 ) ;
switch numel(I) ; 
	case 1
		[dudz] = gradient(u,z) ;
		[dvdz] = gradient(v,z) ;
	case 2
		[~,dudz] = gradient(u,1,z) ;
		[~,dvdz] = gradient(v,1,z) ;
	case 3
		[~,dudz,~] = gradient(u,1,z,1) ;
		[~,dvdz,~] = gradient(v,1,z,1) ;
end

S = sqrt( dudz.^2 + dvdz.^2 ) ;

end
