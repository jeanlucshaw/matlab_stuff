function [ue,F] = currextrap(u,N2,p,z,varargin) 
% Synthax:          [ue,F]  = currextrap(u,N2,p,z,varargin)
% 
% Uses the dynmodes routine and genetic optimisation to fit the first
% 5 vertical modes on current profile 'u', where 'N2' is the buoyancy
% frequency, 'p' is the pressure, and 'z' is the depth (positive) of
% the associated current profile. 
%
% Output 'ue' is the current profile fit provided and 'F' is a 5
% element vector representing the weights of the 4 first modes plus
% the speed offset. This is the solution set of the genetic optimising
% routine.
% 
% Optional arguments should be supplied as parameter/value pairs. 
% Supported options are:
%
%	- section,[string]	Give higher importance to fitting the
%				surface 'top', bottom 'bottom', or the
%				entire velocity profile 'all'. In
%				effect this applies a weighing function
% 				which decreases as zc/z, increases as
% 				z/zc, where zc is a constant,
%				or is flat with respect to depth
%				to the compute the scoring function.
%
%	- zc,[float]		Depth where the weighing function becomes
%				unity. Defaults to 1.
%
%	- itermax,[int]		Set the maximum number of iterations
%				of the genetic algorithm. Defaults to
%				500.
%
%	- popsize,[int]		Set the population size for the genetic
%				algorithm. Defaults to 20.
%
%	- nsol,[int]		Number of solution sets computed and
%				from the genetic algorithm. Only the
%				highest scoring solution is returned
%				Defaults to 1.  
% Dependecies:
%	- dynmodes.m	: https://github.com/sea-mat/dynmodes
% 	- genopt.m	: https://github.com/jeanlucshaw/matlab_stuff/blob/master/math/genopt.m   


%OPTIONS	

%default
zc = 1 ;
NP = 20 ; 	% Genetic algo population size
NI = 500 ; 	% Genetic algo max iterations
NS = 1 ;	% Genetic algo number of solutions averaged

%user
for ii = 1:2:nargin-5
	switch varargin{ii}
		case 'section'
			switch varargin{ii+1}
				case 'top'
					W	= zc./z;
				case 'bottom'
					W	= z./zc;
				case 'all'
					W	= ones(size(z));
			end
		case 'itermax'
			NI = varargin{ii+1} ;
		case 'popsize'
			NP = varargin{ii+1} ;
		case 'nsol'
			NS = varargin{ii+1} ;
		case 'zc'
			zc = varargin{ii+1} ;
	end
end

I  = find(isfinite(N2)) ;
[wmodes,pmodes,ce]=dynmodes(N2(I)',p(I)',1) ;

%REGRID
pm = interp1(1:0.5:1+0.5*numel(I),pmodes,z,'linear','extrap') ;

%FIT
I	= isfinite(u) ;
Q	= @(a) sum( W(I).*sqrt( (a(1)*pm(I,1)+a(2)*pm(I,2)+a(3)*pm(I,3)+a(4)*pm(I,4)+a(5)*pm(I,5)+a(6) - u(I)).^2 ) ) ;
A	= nan([NS,6]) ;
HS	= nan([NS,1]) ;
for ii = 1:NS
	[A(ii,:),HS(ii)]   = genopt(Q,0.01,[-0.2 0.2;-0.2 0.2;-0.2 0.2;-0.2 0.2;-0.2 0.2;-2 2],NP,NI,0.8,0.2) ;
end

%OUTPUT
[~,I]	= min(HS) ;
F	= A(I,:) ;
ue	= F(1)*pm(:,1)+F(2)*pm(:,2)+F(3)*pm(:,3)+F(4)*pm(:,4)+F(5)*pm(:,5)+F(6) ;
end
