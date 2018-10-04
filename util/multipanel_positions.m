function output = multipanel_positions(I,J,bp,tp,rp,lp,vs,hs)
% Synthax : output = multipanel_positions(I,J,bp,tp,rp,lp,vs,hs)
%                    
% Returns a matrix 'output' of size [Npan 4] where 'Npan' is the number of
% pannels for a given figure. Intended use is the feed as input to axes.
%
% Inputs:
%	
%	- I	: number of rows
%	- J	: number of columns
%	- bp	: portion of vs used as padding at the bottom	[0-1]
%	- tp	: portion of vs used as padding at the top	[0-1]
%	- rp	: portion of hs used as padding to the right	[0-1]
%	- rp	: portion of hs used as padding to the left	[0-1]
%	- vs	: total vertical white space			[0-1]
%	- hs	: total horizontal white space			[0-1]
%	
Npan = I*J ;

hp    = (hs-rp-lp)/(J+1) ;	% horizontal pad
hv    = (vs-bp-tp)/(I+1) ;	% vertical   pad

h     = (1-vs)/I  ;		% vertical   panel size
w     = (1-hs)/J ;		% horizontal panel size

output = repmat(nan,[Npan 4]) ;
row    = 0 ;

K = 1;
for jj = 1:J
	for ii=1:I
		output(K,1) = hp + lp + (jj-1)*(w+hp);
		output(K,2) = 1  - (h+hv)*(ii) - tp;
		output(K,3) = w;
		output(K,4) = h;
	K = K + 1;
	end
end 

end
