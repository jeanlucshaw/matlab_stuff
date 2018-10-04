function torplot(theta,phi,varargin)
% Synthax :            torplot(theta,phi,varargin)
%
% Plot two angular variables 'theta' and 'phi' against eachother projected
% on a toroidal space. Optional parameters must be entered as parameter
% value pairs. Available options are:
%
% 	labels,[cell array]	contains the two strings labeling the axes
%
%	offsets,[vector array] 	three value vector defining the angle at
%				which will be draw respectively the
%				starts of the theta axis, the start of the
%				phi axis, and the start of the one to one
%				line. They all default to zero. 
%
%	degrad,[string]		if value is 'deg', inputs will be converted
%				from degrees to radians. Otherwise, inputs
%				will be assumed to be radians.
%	
%	colormap,[nx3 matrix]	colormap to be used on the scatter plot. If no
%				value is provided dots will be red.
%
%	intensity,[px1 vector]	value to be displayed by the color mapping. 'p'
%				is the length of theta and phi. If no value is
%				provided all dots will be red.
%
%	gridon,[string]		display grid in theta direction. If 'y' is provided
%				grid will be displayed.
%
%	grid,[vector]		theta grid values. Defaults to [90 180 270].
%
%	cdir,[string]		direction in which the outside axe will run.
%				To make the ax run clockwise set this parameter
%				to 'cw'. The default is to run counter clockwise.
%				This can be specified by the string 'ccw'.
%

%% OPTIONS
% default
thelab = '' ;
philab = '' ;
thes   = 0 ;
phis   = 0 ;
dots   = 0 ;
degrad = 'rad' ;
phiang = 45;
grid   = [90 180 270] ;
cdir	= 'ccw' ;

% user
for ii = 1:2:numel(varargin)-1 
	switch varargin{ii}
		case 'labels'
			thelab = varargin{ii+1}(1) ;
			philab = varargin{ii+1}(2) ;
		case 'offsets'
			thes   = varargin{ii+1}(1) ;
			phis   = varargin{ii+1}(2) ;
			dots   = varargin{ii+1}(3) ;
		case 'degrad'
			degrad = varargin{ii+1} ;
		case 'colormap'
			cmap   = varargin{ii+1} ;
		case 'intensity'
			int    = varargin{ii+1} ;
		case 'cdir'
			cdir   = varargin{ii+1} ;
	end
end

% CONVERT TO RADIANS
if startsWith(degrad,'deg')
	theta = theta*pi/180 ;
	phi   = phi*pi/180;
end

% MAKE COLOR VECTOR
if exist('cmap') & exist('int')
	cols	= nan([numel(int) 3]) ;
	nmap	= numel(cmap(:,1)); 
	I	= 1 + floor( (nmap-1)*int/max(int) ); 
	for ii  = 1:numel(int)
		cols(ii,:)	= cmap(I(ii),:) ;
	end
end

% DRAW STUFF
figure; 
axes('position',[0 0 1 1])
r0 = 0.2 ;
rM = 1 ;
r  = rM - r0 ;

% CONVERT INPUT TO TOROIDAL COORDINATES
oneone  = [1:360]*pi/180 ;
if strcmp(cdir,'cw')
	% input data
	x  = (r0 + r*theta/(2*pi)).*sin(phi) ;
	y  = (r0 + r*theta/(2*pi)).*cos(phi) ;
	% gid
	gx	= nan(numel(grid),1) ;
	gy	= nan(numel(grid),1) ;
	for ii	= 1:numel(grid)
		gx(ii)	= 1.1*rM*sind(grid(ii)) ;	
		gy(ii)	= 1.1*rM*cosd(grid(ii)) ;	
	end
	% one to one ratio
	x11  = (r0 + r*oneone/(2*pi)).*sin(oneone) ;
	y11  = (r0 + r*oneone/(2*pi)).*cos(oneone) ;
else
	% input data
	x  = (r0 + r*theta/(2*pi)).*cos(phi) ;
	y  = (r0 + r*theta/(2*pi)).*sin(phi) ;
	% grid
	gx	= nan(numel(grid),1) ;
	gy	= nan(numel(grid),1) ;
	for ii	= 1:numel(grid)
		gx(ii)	= 1.1*rM*cosd(grid(ii)) ;	
		gy(ii)	= 1.1*rM*sind(grid(ii)) ;	
	end
	% one to one ratio
	x11  = (r0 + r*oneone/(2*pi)).*cos(oneone) ;
	y11  = (r0 + r*oneone/(2*pi)).*sin(oneone) ;
end

if exist('cmap') & exist('int')
	scatter(x,y,2,cols,'filled') ; hold on;
	colormap(cmap) ;
	caxis([min(int) max(int)])
	jlcolorbar('$u_w$ (m/s)','location','east','position',[0.85 0.2 0.02 0.6]) ;
else
	scatter(x,y,2,'r','facecolor','r') ; hold on;
end

% PLOT TOR AXES
plot(r0*cosd([1:360 1]),r0*sind([1:360 1]),'k') ; 
plot(rM*cosd([1:360 1]),rM*sind([1:360 1]),'k') ; 

% PLOT GRID
rgrid = r0+r*grid/360 ;
gtextang = 45 ;
for ii =1:numel(grid)
	plot(rgrid(ii).*cosd([5:355]+gtextang),rgrid(ii).*sind([5:355]+gtextang), ...
                            'color',[0.5 0.5 0.5], ...
                            'linestyle','--') ; 
	text(rgrid(ii).*cosd(gtextang),rgrid(ii).*sind(gtextang),num2str(grid(ii)), ...
             'rotation',gtextang-90, ... 
             'horizontalalignment','center', ...
             'fontsize',9, ...
             'color',[0.5 0.5 0.5]) ;
end
%text(1.1,0,'0','horizontalalignment','center') ;
%text(0,1.1,'90','horizontalalignment','center') ;
%text(-1.1,0,'180','horizontalalignment','center') ;
%text(0,-1.1,'270','horizontalalignment','center') ;
for ii	= 1:numel(gx)
	text(gx(ii),gy(ii),num2str(grid(ii)),'horizontalalignment','center') ;
end


% PHI
phiax = [0:10]+phis-7.5+phiang ; 
plot(cosd(phiax),sind(phiax),'k','linewidth',3) ;
patch([cosd(phiax(end) + 5),0.95*cosd(phiax(end)),1.05*cosd(phiax(end))], ...
      [sind(phiax(end) + 5),0.95*sind(phiax(end)),1.05*sind(phiax(end))],'k') ;
plot(cosd(phis),sind(phis),'ko','markerfacecolor','k')
text(1.1*cosd(phis+phiang),1.1*sind(phis+phiang),philab,'horizontalalignment','center') ;
%text(1.1*cosd(phis-4),1.1*sind(phis-4),'0','horizontalalignment','center') ;

% THETA 
wid   = 0.05;
theax = [r/2+r0/2  r/2+7*r0/6] ; 
plot(theax*cosd(thes),theax*sind(thes),'k','linewidth',3) ;
patch([(theax(2)+0.09)*cosd(thes),theax(2)*cosd(thes)-wid*sind(thes),theax(2)*cosd(thes)+wid*sind(thes)], ...
      [(theax(2)+0.09)*sind(thes),theax(2)*sind(thes)+wid*cosd(thes),theax(2)*sind(thes)-wid*cosd(thes)],'k') ;
plot([r0*cosd(thes) rM*cosd(thes)],[r0*sind(thes) rM*sind(thes)],'k')
%text(theax(2)*cosd(thes+15),theax(2)*sind(thes+15),thelab,'horizontalalignment','center') ;
plot(r0*cosd(thes),sind(thes),'ko','markerfacecolor','k')
text(0.5*r0*cosd(thes),0.5*r0*sind(thes),thelab,'horizontalalignment','center') ;

% ONE TO ONE LINE
plot(x11,y11,'k')

% GRAPHICAL PARAMETERS
box off;
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'xticklabels',[],'xtick',[],'xcolor','none')
set(gca,'yticklabels',[],'ytick',[],'ycolor','none')
if exist('cmap') & exist('int')
	axis([-1.3 1.7 -1.1 1.1])
else
	axis([-1.3 1.3 -1.3 1.3])
end

end
