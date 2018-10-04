function top	= listmax(x,N,frmt)
% Syntax : 	top	= listmax(x,N,frmt)
%
% Display the 'N' maximum values from the first column
% of 'x' in the string format 'frmt'. If 'x' is a matrix
% with M columns, listmax will sort the whole matrix 
% according to the rows of the first column. 'frmt' needs 
% to contain M spots for string interpolation. The 'N' top
% values of the sorted matrix are returned as 'top'. 'x' may
% contain NaNs but they are omitted.
%
% 

% Make sorted, de-nan-ed list of input
[~,I]	= sort(x(:,1),'descend') ;
x	= x(I,:) ;
I	= find( isnan(x(:,1)) ) ;
x(I,:)	= [] ;

% Print the results
for ii	= 1:N
	disp(sprintf(frmt,x(ii,:))) ;
end

% Output the top
top	= x(1:N,:) ;
end
