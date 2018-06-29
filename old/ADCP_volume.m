function ADCP_volume(ADCP_file_name,theta,corrtag,WD)
disp('========== 3) Time interpolation ==========') ;

CC      = strrep(num2str(theta),'.','') ;

load([WD 'data/ADCP_23237/' ADCP_file_name(1:end-4) '_gridded_' CC '_' corrtag '.mat']) ;

%% INITIALISE OUTPUT STRUCTURE
V.along = [] ;
V.cross = [] ;
V.lon   = C.lon ;
V.lat   = C.lat ;
V.dist  = C.dist ;
V.z     = C.z ;

%% SET TIME
daystr = datestr(min(C.time),'yyyy-mm-dd');
daystrstrt = [daystr 'T11:00'];
daynumstrt = datenum(daystrstrt,'yyyy-mm-ddTHH:MM') ;
daynumstop = daynumstrt + 13./24 ;                        % 13 hours later
tstep   = 10./(24*60*60) ; % 10 seconds (20 pings)
%V.time  = min(C.time):tstep:max(C.time) ; 
V.time  = daynumstrt:tstep:daynumstop ; 


%% INTERPOLATE 
for ii = 1:numel(C.grdist)
    I  = find(V.dist == C.grdist(ii)) ;
    sheet = C.along(:,I) ;
    stime = C.time(I) ;
    isheet = interp1(stime,sheet',V.time)' ;
    
    V.along = cat(3,V.along,isheet) ;
    
    sheet = C.cross(:,I) ;
    isheet = interp1(stime,sheet',V.time)' ;
    
    V.cross = cat(3,V.cross,isheet) ;
end

V.grdist = C.grdist ;
V.suffix = C.suffix ;

% figure;
% 
% TTT = 3150 ;
% 
% [C,h] = contourf(squeeze(V.along(TTT,:,:)),[-0.25:0.01:0.25]) ;
% hold on 
% [C2,h2]=contour(squeeze(V.along(TTT,:,:)),[0 0]) ;
% set(h,'linestyle','none')
% set(h2,'linecolor','k')
% axis ij
% caxis([-0.25 0.25])
% colorbar;

figure('visible','off')
pcolor(V.time,V.grdist,squeeze(mean(V.along,1,'omitnan'))')
shading flat;
hold on
plot(C.time,C.dist,'k-','linewidth',2)
datetick('x','HH:MM')
xlabel('Time (UTC)')
ylabel('x$_T$ (km)')
axis([min(C.time) max(C.time) min(C.dist) max(C.dist)])
colorbar;

paperwidth = 20 ;
paperheight = 12 ;
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/' ... 
             '/ADCP_cleanup/volume_along_TV_da_' V.suffix(2:end) ...
             '_' CC '_' corrtag '.png']);

figure('visible','off')
pcolor(V.time,V.grdist,squeeze(mean(V.cross,1,'omitnan'))')
shading flat;
hold on
plot(C.time,C.dist,'k-','linewidth',2)
datetick('x','HH:MM')
xlabel('Time (UTC)')
ylabel('x$_T$ (km)')
axis([min(C.time) max(C.time) min(C.dist) max(C.dist)])
colorbar;

paperwidth = 20 ;
paperheight = 12 ;
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/' ... 
             '/ADCP_cleanup/volume_cross_TV_da' V.suffix(2:end) ...
             '_' CC '_' corrtag '.png']);

save([WD 'data/ADCP_23237/' ADCP_file_name(1:end-4) '_volume_' CC '_' corrtag ],'V') ;

end
