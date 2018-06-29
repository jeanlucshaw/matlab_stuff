function myprint(path,pw,ph)
% Synthax :               myprint(path,pw,ph)
% 
% Prints current figure with options -dpng and -r300. Saves the file in
% 'path' at paperwidth 'pw' and paperheight 'ph'.
%

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[pw pw]);
set(gcf,'PaperPosition',[0 0 pw ph]);
print('-dpng', '-r300',path);  
end