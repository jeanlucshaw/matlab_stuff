function output = unifrnd(b_min,b_max,d_len,m,n)

% Synthax :          output = unifrnd(b_min,b_max,d_len,m,n)
%
% Imitates the unifrnd function from the statistics toolbox. Generates a
% 'm' by 'n' matrix of values chosen from a uniform distribution of minimum
% 'b_min' and maximum 'b_max' . The parameter 'd_len' controls the number
% of elements in the uniform distribution the 'output' values are chosen
% from.
%

dist   = linspace(b_min,b_max,d_len)' ;
I      = randi([1 d_len],[m n]) ;

output = dist(I) ;

end