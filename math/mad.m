function Y	= mad(X,varargin)
% Syntax : Y	= mad(X,varargin)
%
% Absolute deviation from the median. Runs the following one-liner.
%
% Y	= median( abs(X - median(X,varargin{:})),varargin{:}) ;
%
	Y	= median( abs(X - median(X,varargin{:})),varargin{:}) ;
end	
