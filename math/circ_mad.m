
function Y	= mad(X,varargin)
% Syntax : Y	= mad(X,varargin)
%
% Absolute deviation from the circular median. Runs the following one-liner.
%
% Y	= circ_median( abs(X - circ_median(X,varargin{:})),varargin{:}) ;
%
% This function is complementary to circ_median, for the cstats Matlab package
% written by Phillip Berens.
%
	Y	= circ_median( abs(X - circ_median(X,varargin{:})),varargin{:}) ;
end	
