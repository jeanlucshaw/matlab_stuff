function ADCP_gridding(ADCP_file_name,G,Gstrt,Gstop,Fh,Fv,theta,corrtag,WD)

disp('========== 2) Smoothing then griding ==========')

CC             = strrep(num2str(theta),'.','') ; 

load([WD 'data/ADCP_23237/' ADCP_file_name(1:end-4) '_CC' CC 'd_' corrtag '.mat']) ;

paperwidth  = 20 ;
paperheight = 6 ;

suffix = ['_1DF' num2str(Fh) 'm_G' sprintf('%d',1000*G) 'm'] ;

%% TRANSFORM FILTER SIZE -> NUMBER OF CASTS
hsize = round(Fh./(mean(diff(S.time)*24*3600,'omitnan').*2));

%% INITIALISE output structure
C.along = [] ;
C.cross = [] ;
C.dist  = [] ;
C.time  = [] ;
C.lon   = [] ;
C.lat   = [] ;
sm_al   = [] ;
sm_cr   = [] ;
sm_time = [] ;

%% Find index of extreme sides of transects
[~,MAXs] = findpeaks( S.dist,'MinPeakHeight',3.7,'MinPeakProminence',0.01) ;
[~,MINs] = findpeaks(-S.dist,'MinPeakHeight',-0.3,'MinPeakProminence',1) ;
II       = sort([MAXs MINs]) ;
II       = [1 II] ;

% figure
% plot(S.time,S.dist,'k.')
% hold on
% plot(S.time(II),S.dist(II),'o')
% datetick('x','HH:MM:SS')

%% SET GRID 
grdist    = [Gstrt:G:Gstop]' ;

for ii = 1:numel(II)-1

    start = ii ;
    stop = start + 1 ;

    % % VIEW THIS TRANSECT
    % figure ;
    % pcolor(S.time(II(start):II(stop)),S.z,S.cross(:,II(start):II(stop)))
    % hold on
    % datetick('x','HH:MM')
    % xlabel('Time')
    % ylabel('Depth (m)')
    % shading flat ;
    % axis ij;
    % ylim([0 80]) ; xlim([S.time(II(start)) S.time(II(stop))]) ;
    % caxis([-1 1])
    % jlcolorbar('u$_{cross}$') ;
    
    %% SELECT DATA FOR THIS TRANSECT
    DIST   = S.dist(II(start):II(stop))' ;
    LON    = S.plon(II(start):II(stop)) ;
    LAT    = S.plat(II(start):II(stop)) ;
    TIME   = S.time(II(start):II(stop)) ;
    ALONG  = S.along(:,II(start):II(stop)) ;
    CROSS  = S.cross(:,II(start):II(stop)) ;
    pNAN   = S.pNanUnder(:,II(start):II(stop)) ;
    
    %% MAKE THE DATA REGULAR
    if DIST(1) > DIST(end)
        % Exclude data too close to boat turns
        I = find(DIST <= Gstop + 0.1 & DIST >= Gstrt - 0.1 ) ;
        DIST        = DIST(I) ;
        LON         = LON(I) ;
        LAT         = LAT(I) ;
        TIME        = TIME(I) ;
        ALONG       = ALONG(:,I) ;
        CROSS       = CROSS(:,I) ;
        pNAN        = pNAN(:,I) ;
        
        % Sort according to transect direction
        [DIST,Is]   = sort(DIST,'descend') ;
        % Clip double values
        [DIST,IA,~] = unique(DIST,'stable') ;
        TIME        = TIME(IA) ;
        LON         = LON(IA) ;
        LAT         = LAT(IA) ;
        ALONG       = ALONG(:,IA) ;
        CROSS       = CROSS(:,IA) ;
        pNAN        = pNAN(:,IA) ;
    else
        % Exclude data too close to boat turns
        I = find(DIST <= Gstop + 0.1 & DIST >= Gstrt - 0.1) ;
        DIST        = DIST(I) ;
        TIME        = TIME(I) ;
        LON         = LON(I) ;
        LAT         = LAT(I) ;
        ALONG       = ALONG(:,I) ;
        CROSS       = CROSS(:,I) ;
        pNAN        = pNAN(:,I) ;        

        % Sort according to transect direction
        [DIST,Is]   = sort(DIST,'ascend') ;
        % Clip double values
        [DIST,IA,~] = unique(DIST,'stable') ;
        TIME        = TIME(IA) ;
        LON         = LON(IA) ;
        LAT         = LAT(IA) ;
        ALONG       = ALONG(:,IA) ;
        CROSS       = CROSS(:,IA) ;
        pNAN        = pNAN(:,IA) ;
    end
    
    %% SMOOTH DATA
    I        = find( pNAN >= 0.99) ;
    ALONG    = movmean(ALONG,hsize,2,'omitnan') ;
    ALONG    = movmean(ALONG,Fv,1,'omitnan');
    ALONG(I) = nan ;
    
    I        = find( pNAN >= 0.99) ;    
    CROSS    = movmean(CROSS,hsize,2,'omitnan') ;
    CROSS    = movmean(CROSS,Fv,1,'omitnan') ;
    CROSS(I) = nan ;
    
    % Save the smooth uninterpollated data to inspect
    sm_al   = [sm_al repmat(nan,[numel(ALONG(:,1)) 1]) ALONG(:,2:end)] ;
    sm_cr   = [sm_cr CROSS] ;
    sm_time = [sm_time ; TIME] ;
    
    %% INTERPOLATE ON SPATIAL GRID with TIME TAGS
    if DIST(1) > DIST(end)
        tdist   = flipud(grdist) ;
        grtime  = interp1(DIST,TIME,tdist) ;
        Inonan  = find(~isnan(grtime)) ;
        grtime  = grtime(Inonan) ;
        tdist   = tdist(Inonan) ;
        gralong = interp1(DIST,ALONG',tdist)' ;
        grcross = interp1(DIST,CROSS',tdist)' ;
        grlon   = interp1(DIST,LON,tdist) ;
        grlat   = interp1(DIST,LAT,tdist) ;
    else
        tdist   =        grdist  ;
        grtime  = interp1(DIST,TIME,tdist) ;
        Inonan  = find(~isnan(grtime)) ;
        grtime  = grtime(Inonan) ;
        tdist   = tdist(Inonan) ;
        gralong = interp1(DIST,ALONG',tdist)' ;
        grcross = interp1(DIST,CROSS',tdist)' ;
        grlon   = interp1(DIST,LON,tdist) ;
        grlat   = interp1(DIST,LAT,tdist) ;
    end
    
    
    %% SAVE TO THE OUTPUT STRUCTURE
    C.time  = [C.time ; grtime] ;
    C.dist  = [C.dist ; tdist]  ;
    C.lon   = [C.lon  ; grlon ] ;
    C.lat   = [C.lat  ; grlat ] ;
    C.along = [C.along gralong] ;
    C.cross = [C.cross grcross] ;

end

%% Not sure why i need this but i do
% [C.time,IA,~] = unique(C.time,'stable') ;
% C.dist        = C.dist(IA) ;
% C.lon         = C.lon(IA) ;
% C.lat         = C.lat(IA) ;
% C.along       = C.along(:,IA) ;
% C.cross       = C.cross(:,IA) ;

C.z      = S.z ;
C.grdist = grdist ;
C.suffix = suffix ;

save(['data/ADCP_23237/' ADCP_file_name(1:end-4) '_gridded' '_' CC '_' corrtag],'C') ;

%% PRINT THE SMOOTHED DATA
figure('visible','off');

pcolor(sm_time,C.z,sm_al)
hold on
plot(S.time,S.depth*cosd(25),'k')
plot(S.time,S.depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(S.time) max(S.time)]) ;
caxis([-1 1])
jlcolorbar('u$_{along}$') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/along_smoo' suffix '_' CC '_' corrtag '.png']);

figure('visible','off');

pcolor(sm_time,C.z,sm_cr)
hold on
plot(S.time,S.depth*cosd(25),'k')
plot(S.time,S.depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(S.time) max(S.time)]) ;
caxis([-1 1])
jlcolorbar('u$_{cross}$') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/cross_smoo' suffix '_' CC '_' corrtag '.png']);

% figure ;
% pcolor(C.time,C.z,C.along)
% hold on
% datetick('x','HH:MM')
% xlabel('Time')
% ylabel('Depth (m)')
% shading flat ;
% axis ij;
% ylim([0 30]) ; xlim([C.time(1) C.time(end)]) ;
% caxis([-1 1])
% jlcolorbar('u$_{cross}$') ;

figure('visible','off')
for ii=1:numel(grdist)
    plot([min(S.time) max(S.time)],[grdist(ii) grdist(ii)],'-','color',[0.7 0.7 0.7])
    hold on
end
plot(S.time,S.dist,'k','linewidth',1)
plot(C.time,C.dist,'r','linewidth',1)
plot(S.time(II),S.dist(II),'ko')
ylabel('x$_T$ (km)')
xlabel('Time')
datetick('x','HH:MM')
axis([min(S.time) max(S.time) min(S.dist) max(S.dist)])

paperwidth = 20 ;
paperheight = 12 ;
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/sorting' suffix '_' CC '_' corrtag '.png']);

end
