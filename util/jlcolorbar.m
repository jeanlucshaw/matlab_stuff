function cbar = jlcolorbar(S,varargin) 
%
% Synthax : cbar = jlcolorbar(S,varargin) 
%
% Puts a colorbar with ticks and label 'S' read by the  latex interpreter
% on the current axes and returns its handle 'cbar'. I can't believe I
% actually need to write this...

cbar = colorbar('ticklabelinterpreter','latex',varargin{:});
ylabel(cbar,S,'interpreter','latex','fontsize',10)



end
