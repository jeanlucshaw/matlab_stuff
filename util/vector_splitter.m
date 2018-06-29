function cell_array = vector_splitter(input,delimiter)
% Synthax : cell_array = vector_splitter(input,delimiter)
%
% Returns a cell array built from the input vector cut
% at every occasion of the delimiter.

if isnan(delimiter)
    I = find(isnan(input)) ; 
else
    I = find(input==delimiter) ; 
end

cell_array = [] ;

% Premier vecteur
cell_array = [cell_array {input(1:I(1)-1)'}] ;
% Vecteurs centraux
for i=1:numel(I)-1
    cell_array = [cell_array {input(I(i)+1:I(i+1)-1)'}] ;
end
% Dernier vecteur
cell_array = [cell_array {input(I(end)+1:end)'}] ;
end
