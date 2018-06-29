function output = multipanel_positions(I,J,bottompad,rightpad,varargin)
% Synthax : output = multipanel_positions(I,J,bottompad,rightpad,varargin)
%                    
% Returns a matrix 'output' of size [Npan 4] where 'Npan' is the number of
% pannels for a given figure. Intended use is the feed as input to axes. 25
% percent of the figure is used for spacing between the figure sides and
% the pannels and spacing between the pannels themselves. A percentage
% 'bottompad' of the vertical space may be added to the bottom of the
% figure.

Npan = I*J ;

hspace = 0.25.*varargin{1} ;
vspace = 0.25.*varargin{2} ;

hp    = hspace/(J+1) ; % horizontal pad
hv    = vspace/(I+1) ; % vertical   pad

h     = (1-vspace-bottompad)/I  ;     % vertical   panel size
w     = (1-hspace-rightpad) /J ;      % horizontal panel size

output = repmat(nan,[Npan 4]) ;
row    = 0 ;

K = 1;
for jj = 1:J
	for ii=1:I
		output(K,1) = hp + (jj-1)*(w+hp);
		output(K,2) = 1  - (h+hv)*(ii) ;
		output(K,3) = w;
		output(K,4) = h;
	K = K + 1;
	end
end 


%if Npan~=5
%    for kk = 1:Npan
%        output(kk,1) = hp + (kk-1)*(w+hp) ;
%        output(kk,2) = 1 - (hv + h + row*(h+hv)) ;
%        output(kk,3) = w ;
%        output(kk,4) = h ;
%        if mod(kk,J) == 0 ; 
%            row = row + 1 ;
%        end
%    end
%else
%    for kk = 1:4
%        output(kk,1) = hp + (mod(kk+1,J))*(w+hp) ;
%        output(kk,2) = 1 - (hv + h + row*(h+hv)) ;
%        output(kk,3) = w ;
%        output(kk,4) = h ;
%        if mod(kk,J) == 0 ; 
%            row = row + 1 ;
%        end
%    end
%        output(5,1) = (1-w)/2 ;
%        output(5,2) = 1 - (hv + h + row*(h+hv)) ;
%        output(5,3) = w ;
%        output(5,4) = h ;
%end
end
