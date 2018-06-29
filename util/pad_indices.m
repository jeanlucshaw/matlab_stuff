function output = pad_indices(input,psize,maxN)
%
% Synthax : output = pad_indices(input,psize,maxN)
%
% Pads the values of an index vector 'input' with 
% adjacent indices. Pads of size 'psize' are inserted
% on both sides. 'maxN' is the size of the vector from
% which the indices are extracted.

matp = repmat(0,[numel(input) psize]) ;
matn = repmat(0,[numel(input) psize]) ;
for i = 1:psize
    matp(:,i) = input + i ;
end
for i = 1:psize
    matn(:,i) = input - i ;
end

output = unique([matn(:) ; input ; matp(:)]) ;

%CLEAR NEGATIVE IND
ind         = find(output<1) ;
output(ind) = [] ;

%CLEAR INDICES GREATER THAN NUMEL(INPUT)
ind         = find(output>maxN) ;
output(ind) = [] ;

end
