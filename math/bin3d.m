function [vals,sds]	= bin3d(in,bsz,varargin)
% Syntax :	[vals,sds]	= bin3d(in,bsz,varargin)
% 
% Bin the 3D matrix in. By default 'bsz' is a three element
% vector containing the number of bins in each dimension. By
% specifying the parameter value pair 'bintype','npts' this 
% can be changed to the number of elements per bin in each
% dimension. The averaged values are returned in 'vals' and
% their standard deviations are returned in 'sds'.
%

% SIZE OF INPUT
[M,N,L]	= size(in) ;

% DEFAULT OPTIONS
osz	= bsz ;	% Default bintype is number of divisions
bsx	= floor(M/bsz(1)) ;
bsy	= floor(N/bsz(2)) ;
bsz	= floor(L/bsz(3)) ;

% USER OPTIONS
if ~isempty(varargin{:})
for ii	= 1:2:numel(varargin{:})
	switch varargin{ii}
		case 'bintype'
			switch varargin{ii+1}
				case 'ndiv'
					bsx	= floor(M/bsz(1)) ;
					bsy	= floor(N/bsz(2)) ;
					bsz	= floor(L/bsz(3)) ;
					osz	= bsz ;
									
				case 'npts'
					bsx	= bsz(1) ;
					bsy	= bsz(2) ;
					bsz	= bsz(3) ;
					osz	= [floor(M/bsz(1)),floor(N/bsz),floor(L/bsz)] ;
			end
	end
end
end

% INIT OUTPUT
vals	= nan(osz) ;
sds	= nan(osz) ;

% CONDUCT AVERAGING
for ii	= 1:osz(1)
	for jj	= 1:osz(2)
		for kk	= 1:osz(3)
			sx		= 1+(ii-1)*bsx ;
			sy		= 1+(jj-1)*bsy ;
			sz		= 1+(kk-1)*bsz ;
			binvals		= in(sx:ii*bsx,sy:jj*bsy,sz:kk*bsz) ;
			vals(ii,jj,kk)	= mean(binvals(:),'omitnan') ;
			sds(ii,jj,kk)	=  std(binvals(:),'omitnan') ;
		end
	end
end


end
