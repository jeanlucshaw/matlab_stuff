clear all;

M	= 10 ;
N	= 100;
L	= 5;

A	= repmat([1:M]',[1 N L]) ;
B	= repmat([1:N] ,[M 1 L]) ;
Cp	= nan([1,1,L]) ;
Cp(:)	= 1:L ;
C	= repmat(Cp ,[M N 1]) ;
Zp	= [0:M]' ;
%Zp	= [1:M].^2 ; Zp	= Zp' ;
Z	= repmat(Zp,[1 N L]) ;

[Ax,Ay,Az]	= gradient(A,1,Z,1) ;
[Bx,By,Bz]	= gradient(B,1,Z,1) ;
[Cx,Cy,Cz]	= gradient(C,1,Z,1) ;

u	= repmat(sqrt(2).*[1:M]',[1 N L]) ;
v	= repmat(sqrt(2).*[1:M]',[1 N L]) ;

S	= uvz2S(u,v,Z) ;

function S = uvz2S(u,v,z)
% Synthax:            S = uvz2S(u,v,z)
% 
% Wrapper for the gradient function which takes as input
% horizontal velocities 'u' and 'v' as well as their 
% corresponding depths 'z'. Note that 'u','v', and 'z'
% are all assumed to be of same size. The matrices of u and
% v MUST have the z axis along the first dimension.
%
% Returns the horizontal shear:
%
%      S = sqrt( (du/dz)^2 + (dv/dz)^2 )
%
%D = ndims(u) ;
%M = nan([1 D]) ;
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
