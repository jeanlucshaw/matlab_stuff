function ADCP_cleanup(ADCP_file_path,ADCP_file_name,BATH_file,CTD_file,GPS_file,strtsamp,stopsamp,theta,corrtag,CTDtag,WD,zerolon,zerolat)

disp('========== 1) ADCP_cleanup ==========') ;

tic
%% PARAMETERS
paperwidth  = 20 ;
paperheight = 6 ;

%% OUTPUT FILE TAGS
CC      = ['CC' strrep(num2str(theta),'.','') 'd'] ;

%% LOAD data
disp('Loading data ...')
load([WD ADCP_file_path ADCP_file_name]) ;
bathystruc = load(BATH_file) ;
bathyfield = fieldnames(bathystruc) ;
BATHY      = getfield(bathystruc,bathyfield{:}) ;

% dataset decimation (for testing)
loaded = 1:1:numel(Sens.dnum) ;

%% CREATE GRAPHICS AND DATA FOLDERS
if 7 == exist([WD 'figs/ADCP_cleanup'])
    disp('Graphics folder found.') ;
else
    disp('Creating graphics folder : figs/ADCP_cleanup')
    mkdir([WD 'figs/ADCP_cleanup']);
end 

if 7 == exist([WD 'data'])
    disp('data folder found.') ;
else
    disp('Creating data folder : data')
    mkdir([WD 'data']);
end 

%% temps ADCP
I         = find(Sens.dnum(loaded) < strtsamp | Sens.dnum(loaded) > stopsamp) ;
loaded(I) = [] ;
time      = Sens.dnum(loaded) ;

%% matrices des vitesses

u_raw = Wt.vel(loaded,:,1)' ; % m/s
v_raw = Wt.vel(loaded,:,2)' ; % m/s
w_raw = Wt.vel(loaded,:,3)' ; % m/s
disp('Done')

%% vitesse du bateau
disp('Making boat speed matrices ...')
fid = fopen([WD GPS_file]) ;
    frmt = ['%s%s%f%f%s%s' repmat('%s',[1 23])] ;
    data = textscan(fid,frmt,'headerlines',49,'delimiter',',') ;
fclose(fid) ;

[boat_tim,IA,~] = unique(datenum(cell2mat(data{6}),'yyyy-mm-ddTHH:MM:SS'));
boat_lon = data{4}(IA) ;
boat_lat = data{3}(IA) ;

% calculer la vitesse du bateau
[boat_u,boat_v,boat_spd,boat_head] = llt2spd(boat_lon,boat_lat,boat_tim) ;

boat_acc = sqrt(gradient(boat_u,boat_tim*(24*3600)).^2 + gradient(boat_v,boat_tim*(24*3600)).^2) ;

% interpoler sur le temps de l ADCP
boat_lon  = interp1(boat_tim,boat_lon,time) ;
boat_lat  = interp1(boat_tim,boat_lat,time) ;
boat_u    = interp1(boat_tim,boat_u,time) ;
boat_v    = interp1(boat_tim,boat_v,time) ;
boat_head = interp1(boat_tim,boat_head,time) ;
boat_acc  = interp1(boat_tim,boat_acc,time) ;
boat_spd  = interp1(boat_tim,boat_spd,time) ;

% creer les matrices de vitesse de bateau de la taille de vitesses raw
boat_u_mat = repmat(boat_u',[numel(u_raw(:,1)) 1]) ;
boat_v_mat = repmat(boat_v',[numel(v_raw(:,1)) 1]) ;

%% positions du centre des bins
% une profondeur du poisson de 1 m est supposee.
z  = 1 + Info.cell1 + 0.5*Info.cell ...
     + [0 : Info.cell : 93*Info.cell] ; % ecrit ainsi par souci de clarte
 
%% raw data depths
depth = BATHY(boat_lon,boat_lat) ;
clear BATHY ;

%% PRINT THE BOAT SPEED MATRICES
disp('Printing boat speed matrices ... ')
figure('visible','off');

pcolor(time,z,boat_u_mat)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/u_boat_speed_' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,boat_v_mat)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('v') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/v_boat_speed_' CC '_' corrtag '.png']);
close;
disp('Done.') ;

%% PRINT THE RAW DATA
disp('Printing the raw data ... ')
figure('visible','off');

pcolor(time,z,u_raw)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/u_raw_' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,v_raw)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('v') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/v_raw_' CC '_' corrtag '.png']);
close;
disp('Done.')

%% rotation corrective des vitesses mesurees
% une rotation de 36 degree correspond a deux fois la declinaison
% magnetique. L hypothese est que la valeur de -18.x a ete utilisee au lieu
% de la valeur de  18.x . Suite utilisation d'un script qui test les angles
% de 0 a 360 degres et minimise la norme de u et v, la valeur de rotation
% optimal semble etre 32 degres. Un script qui cherche la difference
% moyenne entre le heading trouve par les vitesses du bateau et celles de l
% adcp produit un angle optimal de 27.3 degres pour sa part. Toutefois, cet
% angle de rotation produit encore des artefacts de vitesse de bateau.
% Tatonnement fournit 25 degres comme l'angle de rotation optimal
disp('Compass correct rotation ...')
u_rot = u_raw.*cosd(theta) - v_raw.*sind(theta) ;
v_rot = u_raw.*sind(theta) + v_raw.*cosd(theta) ;
disp('Done.')

clear u_raw ; clear v_raw ;

%% PRINT THE ROTATED DATA
disp('Printing rotated data ... ')
figure('visible','off');

pcolor(time,z,u_rot)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/u_rot_' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,v_rot)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('v') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/v_rot' CC '_' corrtag '.png']);
close;
disp('Done.')

%% retrait des donnees corr < 64

switch corrtag
    case 'CF64'
        disp('Removing uncorrelated ...')
        I        = find(Wt.corr(loaded,:,1)' < 64 | Wt.corr(loaded,:,2)' < 64 ... 
                      | Wt.corr(loaded,:,3)' < 64 | Wt.corr(loaded,:,4)' < 64 ) ;
        u_rot(I) = nan ;
        v_rot(I) = nan ;
        w_raw(I) = nan ;
        disp('Done.')

        clear Wt ; clear Info ; clear Sens ;

%% PRINT THE DATA AFTER REMOVAL OF UNCORELLATED DATA
        disp('Printing the correlated data ... ')
        figure('visible','off');

        pcolor(time,z,u_rot)
        hold on
        plot(time,depth*cosd(25),'k')
        plot(time,depth,'k--')
        datetick('x','HH:MM')
        xlabel('Time')
        ylabel('Depth (m)')
        shading flat ;
        axis ij;
        ylim([0 80]) ; xlim([min(time) max(time)]) ;
        caxis([-1 1])
        jlcolorbar('u') ;

        set(gcf,'PaperUnits','centimeters');
        set(gcf,'PaperSize',[paperwidth paperheight]);
        set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
        print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/u_corr_' CC '_' corrtag '.png']);
        close;
        
        figure('visible','off');

        pcolor(time,z,v_rot)
        hold on
        plot(time,depth*cosd(25),'k')
        plot(time,depth,'k--')
        datetick('x','HH:MM')
        xlabel('Time')
        ylabel('Depth (m)')
        shading flat ;
        axis ij;
        ylim([0 80]) ; xlim([min(time) max(time)]) ;
        caxis([-1 1])
        jlcolorbar('v') ;

        set(gcf,'PaperUnits','centimeters');
        set(gcf,'PaperSize',[paperwidth paperheight]);
        set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
        print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/v_corr_' CC '_' corrtag '.png']);
        close;
        disp('Done.')
    case 'NCF'
end

%% correction des vitesses pour la vitesse du bateau
disp('Boat speed correction ...')
u = u_rot + boat_u_mat ;
v = v_rot + boat_v_mat ;
w = w_raw ;
disp('Done.')

%% PRINT THE BOAT-CORRECTED DATA
disp('Printing the boat speed corrected data ... ')
figure('visible','off');

pcolor(time,z,u)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/u_boat_' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,v)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('v') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/v_boat_' CC '_' corrtag '.png']);
close;
disp('Done.')


%% REMOVE ADCP DATA FROM Acceleration

switch CTDtag
	case 'CTDout'
fid = fopen([WD CTD_file]);
	frmt = '%s%f%f%f%f' ;
	data = textscan(fid,frmt,'headerlines',1,'delimiter',';');
fclose(fid);

statname = data{1} ;
statstrt = data{4} ;
statstop = data{5} ;

for ii = 1:numel(statname)
	I = find(time < statstop(ii) & time > statstrt(ii)) ;
	u(:,I) = nan ;
	v(:,I) = nan ;
	w(:,I) = nan ;
end
end

NN1 = sqrt(u.^2+v.^2) ;
NN2 = sqrt(u_rot.^2+v_rot.^2) ;
% and speed which are bigger since adding boat speed
I        = find(NN1 > NN2) ;
%u(I) = nan ;
%v(I) = nan ;
%w(I) = nan ;
disp('Done.')

%% FILTER BOAT TURNS

N = 250 ;
B = repmat(1./N,[1 N]) ;
A = 1 ;
boat_head_sm = filtfilt(B,A,boat_head) ;

[~,I1]            = SEtrunc(gradient(boat_head_sm),1) ;
boat_head_clr     = boat_head ;
boat_head_clr(I1) = nan;
[~,I2]            = SEtrunc(gradient(boat_head_clr),0.8) ;
I                 = sort([I1 ; I2]);

%u(:,I) = nan ;
%v(:,I) = nan ;
%w(:,I) = nan ;

%% PRINT THE CTD REMOVED DATA
disp('Printing the CTD removed data ... ')
figure('visible','off');

pcolor(time,z,u)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/u_ctdout_' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,v)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('v') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/v_ctdout_' CC '_' corrtag '.png']);
close;
disp('Done.')

clear u_rot ; clear v_rot ; clear w_raw ; 
clear boat_u_mat ; 
clear boat_v_mat;
clear boat_u ;
clear boat_v ;

%% Build mean transect
disp('Projection along mean transect ...')
rlon = boat_lon ;
rlat = boat_lat ;

mlon = linspace(min(rlon),max(rlon),300) ;
P    = polyfit(rlon,rlat,1) ;
mlat = P(1)*mlon + P(2) ;

% figure
% plot(rlon,rlat,'r.','markersize',0.1)
% hold on
% plot(mlon,mlat,'k.')

%% Projection of the coordinates on the mean transect
a = -P(1) ; b = 1 ; c = -P(2) ;
plon = (b*( b*rlon-a*rlat) - a*c ) ./ ( a*a + b*b ) ;
plat = (a*(-b*rlon+a*rlat) - b*c ) ./ ( a*a + b*b ) ;
disp('Done.')

%% PRINT THE FITTED TRANSECT
disp('Disp printing the fitted transect')
figure('visible','off')
plot(rlon,rlat,'-','color',[0.7 0.7 0.7])
hold on
plot(plon,plat,'k-')
xlabel('Longitude')
ylabel('Latitude')

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/mean_transect' CC '_' corrtag '.png']);
close;
disp('Done.') ;

%% Remove data points further than 200 m from the transect 
disp('Remove points too far from mean transect ... ')
for ii = 1:numel(rlon) 
    proj_dist = m_lldist([rlon(ii) plon(ii)],[rlat(ii) plat(ii)]) ;
end

Igone = find(proj_dist > 0.2) ;
u(Igone,:) = nan ;
v(Igone,:) = nan ;
w(Igone,:) = nan ;
disp('Done.')

%% PRINT DATA AFTER REMOVING POINTS TOO FAR FROM MEAN TRANSECT
disp('Printing the data after removal of casts more than 200m from the mean transect ...') ;
figure('visible','off');

pcolor(time,z,u)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/u_toofar_' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,v)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('v') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/v_toofar_' CC '_' corrtag '.png']);
close;
disp('Done.') ;

%% Rejection from side lobe interference
disp('Removing side-lobe contamination ...')
for ii = 1:numel(plon)
    Igone = find(z >= depth(ii)*cosd(25)) ;
    u(Igone,ii) = nan ;
    v(Igone,ii) = nan ;
    w(Igone,ii) = nan ;
end
disp('Done.')

%% PRINT THE SIDE LOBE CORRECTED DATA
disp('Printing the side lobe corrected data ... ') ;
figure('visible','off');

pcolor(time,z,u)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/u_sidelobe' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,v)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('v') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/v_sidelobe' CC '_' corrtag '.png']);
close;
disp('Done.')

%% Rejection of outliers 
disp('Removing outliers |u , v| > mean + 3 sigma ...')

uu = u.^2 ;
vv = v.^2 ;
NN = sqrt( uu + vv ) ;

figure('visible','off')
histogram(NN(:),[0:0.05:3],'normalization','probability') ;
xlabel('$\sqrt{u^2 + v^2}$ (m/s)')
ylabel('pdf($\sqrt{u^2 + v^2}$)')
xlim([0 3])
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 12]);
set(gcf,'PaperPosition',[0 0 20 12]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/NN_hist' CC '_' corrtag '.png']);
close;

u = SEtrunc(u,3) ;     
v = SEtrunc(v,3) ; 

uu = u.^2 ;
vv = v.^2 ;
NN = sqrt( uu + vv ) ;

figure('visible','off')
histogram(NN(:),[0:0.05:3],'normalization','probability') ;
xlabel('$\sqrt{u^2 + v^2}$ (m/s)')
ylabel('pdf($\sqrt{u^2 + v^2}$)')
xlim([0 3])
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 12]);
set(gcf,'PaperPosition',[0 0 20 12]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/NN_hist_clip' CC '_' corrtag '.png']);
close;
clear uu; clear vv; clear NN;
disp('Done.')


%% PRINT THE OUTLIERS CORRECTED DATA
disp('Printing the data after removal of outliers ... ') ;
figure('visible','off');

pcolor(time,z,u)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/u_outliers_' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,v)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('v') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/v_outliers_' CC '_' corrtag '.png']);
disp('Done.')
close;

%% projection des vitesses // et perp au transect
disp('Projecting to frame of reference ...')
transect_angle = max(mlon) - min(mlon) + 1i*(max(mlat) - min(mlat)) ;
transect_angle = 360*angle(transect_angle)./(2*pi) ;

phi = 90 - transect_angle ;

% % make sure that the heading makes sens
% figure
% plot(boat_lon,boat_lat) ;
% hold on
% plot(xfit,fit,':') ;
% quiver(mean(xfit),mean(fit),0.02*cosd(heading),0.02*sind(heading)) ;

% perpendicular points downstream
% Matrice de rotation inverse pour le changement de referentiel
% ange negatif parce que ref' est tourne en sens horaire
c_perp =  u.*cosd(-phi) + v.*sind(-phi) ;
c_para = -u.*sind(-phi) + v.*cosd(-phi) ;
disp('Done.')

clear u ; clear v ;

% calculate heading according to surface measured velocity
% adcp_head  = range_pm180_2_360(180+atan2d( ... 
%     mean(u_rot(:,:),1,'omitnan'), ...
%     mean(v_rot(:,:),1,'omitnan'))) ;

%% PRINT THE PROJECTED DATA
disp('Printing the projected data ... ')
figure('visible','off');

pcolor(time,z,c_perp)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u$_{along}$') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/along_' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,c_para)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u$_{cross}$') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/cross_' CC '_' corrtag '.png']);
disp('Done.')
close;

%% Create transect coordinates
disp('Creating transect coordinates ...')
%zerolon = min(plon) ;
%zerolat = min(plat) ;

for ii = 1:numel(plon) 
    dist(ii) = m_lldist([zerolon plon(ii)],[zerolat plat(ii)]) ;
end
disp('Done.')

%% CLEAR DATA BELOW 35m DEPTH (Not enough for variability analysis)
disp('Clearing data below 31.84 m ... ') ;
switch corrtag
    case 'CF64'
        c_perp(60:end,:) = [] ;
        c_para(60:end,:) = [] ;
        z(60:end)        = [] ;
    case 'NCF'
end
disp('Done.')

%% PRINT THE PROJECTED TOP 35 m
disp('Printing data after clearing data below 31.84 m ... ')
figure('visible','off');

pcolor(time,z,c_perp)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u$_{along}$') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/along_30m_' CC '_' corrtag '.png']);
close;

figure('visible','off');

pcolor(time,z,c_para)
hold on
plot(time,depth*cosd(25),'k')
plot(time,depth,'k--')
datetick('x','HH:MM')
xlabel('Time')
ylabel('Depth (m)')
shading flat ;
axis ij;
ylim([0 80]) ; xlim([min(time) max(time)]) ;
caxis([-1 1])
jlcolorbar('u$_{cross}$') ;

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[0 0 paperwidth paperheight]);
print('-dpng', '-r300',[WD 'figs/ADCP_cleanup/cross_30m_' CC '_' corrtag '.png']);
close;

%% percentage of nan under
disp('pNanUnder matrix ...')

pNanUnder = repmat(0,size(c_perp)) ; 
for ii = 1:numel(time)
    for jj = 1:numel(z)
        pNanUnder(jj,ii) = numel(find(isnan(c_perp(jj:end,ii))))./numel(c_perp(jj:end,ii)) ;
    end
end
disp('Done.')


%% creation de la structure clean
disp('Saving structure ...')
S.name   = 'Mouth of the BSI transect (12H)' ;
% S.u      = u ;
% S.v      = v ;
S.w      = w ;
% S.cross  = cross_smoo ;
% S.along  = along_smoo ;
S.cross  = c_para ;
S.along  = c_perp ;
S.pNanUnder = pNanUnder ;
%S.angle  = transect_angle ;
% S.lon    = boat_lon ;
% S.lat    = boat_lat ;
S.plon   = plon ;
S.plat   = plat ;
S.dist   = dist ;
S.time   = time ;
S.z      = z ;
S.depth  = depth ;

save([WD 'data/ADCP_23237/' ADCP_file_name(1:end-4) '_' CC '_' corrtag ],'S') ;
disp('All done.')
toc

end
