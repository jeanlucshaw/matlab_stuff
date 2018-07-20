function [yfit,A] = nlinfit(rhs,x,y,pars,varargin)
% 
% Synthax :  [yfit,A] = nlinfit(rhs,x,y,pars,cst1,cst2,...,cstn)
%
% An idea found online and remodeled to make it more modular about how to
% fminsearch.m to perform non-linear, multiple parameter function fitting.
% Variables 'x' and 'y' are the input data, i.e., y = y(x) and must be
% vectors. Variable 'rhs' is a function handle refering to te mathematical
% form of the fitted function. For example, to fit a 1/6 power law, define:
%
%                 rhs = @(a,x) a(1)*x.^(1/6) + a(2) ;
%
% where 'a' is a vector formed of the parameters to fit and 'x' is the
% input independent variable. If variable 'pars' is entered as an integer
% it is interpreted as the number of parameters to fit, and sets the first
% guess for all parameters to zero. If 'pars' is entered as a vector, then
% it is interpreted as the first guess for all parameters.
%
% The function returns the fitted values at 'x' as well as the optimized
% parameter vector 'A'.

% Is pars first guess for parameters or number of parameters
sz = size(pars) ;
if     sz(1) == 1 & sz(2) == 1   
            start = zeros([pars 1]) ;
elseif sz(1) == 1 & sz(2) >  1 | sz(2) == 1 & sz(1) >  1
            start = pars ;
end


% Minimize parameters
         A     = fminsearch(@(a)sum( (y - rhs(a,x,varargin{:})).^2 ),start) ;
         yfit  = rhs(A,x,varargin{:}) ;
end
