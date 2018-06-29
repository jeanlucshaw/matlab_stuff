function OGSL_data_cleanup(fname)

% This function cleans data files from https://ogsl.ca
%
% It performs the following modifications :
%
%          1) Removes all " caracters
%          2) Replaces 't,c'    with 't c'
%          3) Replaces '\d, \d' with '\d \d'
%
% It then proceeds to overwrite the original file
%
tic

fid  = fopen(fname,'r') ;
    header = fgetl(fid) ;
    text   = textscan(fid,'%s','headerlines',0,'delimiter','\n') ;
fclose(fid) ;

text  = regexprep(text{:},'(\d{1}),\s(\d{1})','$1 $2') ; % 3)

for i = 1:numel(text)
    text{i} = strrep(text{i},'"','') ;      % 1)
    text{i} = strrep(text{i},'t,c','t c') ; % 2)
end

% Overwriting the original file
fid = fopen(fname,'w') ;
    fprintf(fid,'%s\n',header ) ;
    fprintf(fid,'%s\n',text{:}) ;
fclose(fid) ;

toc

end