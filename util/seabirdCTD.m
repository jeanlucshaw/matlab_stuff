function varargout = seabirdCTD(path,varargin)
%    Synthax :    varargout = seabirdCTD(path,varargin)
%
% Read Seabird electronics CTD .cnv file(s) automatically, select downcast, 
% and bin average the cast profiles. The data is organised into output 
% structure 'S' with time as field 'dnum' and depth as field 'z'. Other 
% fields depend on the data contained in the raw .cnv file. Detected data
% types are:
% 	
% 	DATA 		field name
%
%	temperature		->	S.temp
%	salinity		->	S.sal
%	fluorescence		->	S.fluo
%	turbidity		->	S.turb
%	oxygen			->	S.oxy
%	longitude		->	S.lon
%	latitude		->	S.lat
%	flag			->	S.flag
%
% Optional parameters should be supplied in parameter/value pairs. Available
% options are:
%
% 	warmup, [integer] 		number of minutes allowed for probe warmup
%					when sampling. Defaults to 0.
%	verticalVelocity, [float]	minimum vertical velocity kept.
%					defaults to 0.
%
% 	binSize, [float]		averaging bin size in (m). 
%					defaults to 1 m.	 
%
%	showz, [string]			plot raw depth against time and cleaned
%					binned depth agains time. 
%					Options are: 'y' or 'n'
%					Defaults to 'n'.
%
% 	namest, [cell array] 		extract information from the file name.
%					The first element of the string array
%					is the delimiter, the remaining are pairs
%					of string and integer
%
%	gpsgar, [string]		path to a Garmin .csv gps track log. Track
%					log is interpolated to the sampling time
%					vector. Probably not general to all Garmin
%					devices. If problematic edit lines 83-92.
%
%	gsw, [string]			calculate absolute salinity, conservative
%					temperature and density using the Gibbs
%					Sea Water matlab package. The package must
%					be on the matlab path for this to work.
%					Added fields to the structure are
%
%					cons. temperature	-> 	S.ctemp
%					cons. density		-> 	S.pden
%					absolute salinity	->	S.asal
%
% If the 'path' variable points to a file, only this file is processed and the
% structure 'S' is returned as a variable. If it points to a folder, all the
% .cnv files present in this folder are processed and the output suppressed.
%


% MANAGE OPTIONS
% defaults
warmup  =  0;
dzthres =  0;
bsz     =  1; 
showz   = 'n' ;
gsw	= false;

% user
if numel(varargin) > 1 
	for ii = 1:2:numel(varargin)
		switch varargin{ii}
			case 'warmup'
				warmup	 = varargin{ii+1}*60 ;
			case 'verticalVelocity'
				dzthres	 = varargin{ii+1} ;
			case 'binSize'
				bsz	 = varargin{ii+1} ;
			case 'showz'
				showz	 = varargin{ii+1} ;
			case 'namestr'	
			Nin = numel(varargin{ii+1}) ;
				delim    = varargin{ii+1}{1}; 
				field 	 = cell([(Nin-1)/2 1]) ;
				namef 	 = cell([(Nin-1)/2 1]) ;
				for jj = 1:(Nin-1)/2
					field{jj}	 = varargin{ii+1}{2*jj} ; 
					namef{jj}	 = varargin{ii+1}{2*jj+1} ; 
				end
			case 'gpsgar' % read navigation from a garmin gps .csv
				fid = fopen(varargin{ii+1}) ;
					ln	 = fgetl(fid);
					while ~startsWith(ln,'ID,trkseg'); ln = fgetl(fid) ;  end
					frmt	 = ['%f%f%f%f%f' repmat('%s',[1 numel(split(ln,','))-5])] ;
					A	 = textscan(fid,frmt,'delimiter',','); 
					lon	 = A{4};
					lat	 = A{3};
					tim	 = datenum(A{6},'yyyy-mm-ddTHH:MM:SS'); 
				fclose(fid);
			case 'gsw'
				if startsWith(varargin{ii+1},'y'); gsw = true; end
		end
	end
end

% FIND + ORGANISE FILE CANDIDATES
if isdir(path)
	dirst     = dir([path '/*.cnv']) ;
	fname     = {dirst.name}';
	fpath     = {dirst.folder}';
else
	dirst     = dir([path]) ;
	fname     = {dirst.name}';
	fpath     = {dirst.folder}';
end

for kk = 1:numel(fname)

% DISPLAY FILENAME
disp(['Processing file : ' fname{kk}])

fid = fopen([fpath{kk} '/' fname{kk}]) ;
ln  = fgetl(fid) ;

sens    = {'temperature','salinity','Depth','oxygen','fluorescence','turbidity','flag','latitude','longitude'} ;
Nsens   = 0 ;
colname = [] ;
while ~startsWith(ln,'*END*')
	ln = fgetl(fid) ;
	% get time related values 
	if startsWith(ln,'# nvalues') ;
		spln    = split(ln); 
		nvalues = str2num(spln{4}); end
	if startsWith(ln,'# interval') ;
		spln    = split(ln) ;
		sr      = str2num(spln{end}); end
	if startsWith(ln,'# start_time') ;
		spln    = split(ln) ;
		st      = datenum([spln{4:7}],'mmmddyyyyHH:MM:SS'); end

	% figure out probe channels
	if startsWith(ln,'# name') ;
		Nsens = Nsens + 1 ;
		tmp   = regexpi(ln,sens,'match') ;
		I     = find(~cellfun('isempty',tmp)) ;
		switch I 
			case 1
				colname = [colname; {'temp'}] ;   
			case 2
				colname = [colname; {'sal'}] ;   
			case 3
				colname = [colname; {'z'}] ;  Iz = Nsens ; 
			case 4
				colname = [colname; {'oxy'}] ;   
			case 5
				colname = [colname; {'fluo'}] ;   
			case 6
				colname = [colname; {'turb'}] ;   
			case 7
				colname = [colname; {'flag'}] ;   
			case 8
				colname = [colname; {'lat'}] ;   
			case 9
				colname = [colname; {'lon'}] ;   
		end
	end	
end

%% READ DATA
A = textscan(fid,repmat('%f',[1 Nsens])) ;

%% MAKE TIME VECTOR 
dnum = st+sr*[0:nvalues-1]'/(24*3600) ;

%% GET DOWNCAST
z      = movmean(A{:,Iz},10/sr) ;
if startsWith(showz,'y'); figure; plot(dnum,z,'b'); hold on; end

dz     = gradient(z,dnum*24*3600) ;					% compute downward velocity 
z( [1:nvalues]' < warmup/sr | dz<=0.001  ) = nan ;			% get rid of start and end values
[~,Im] = min(z) ;							% find min depth 
[~,IM] = max(z) ;							% find max depth 
I      = find(dz > dzthres & [1:nvalues]' > Im & [1:nvalues]' < IM)  ;	% select down-going between min z and max z
if startsWith(showz,'y'); plot(dnum(Im),z(Im),'*'); plot(dnum(IM),z(IM),'o'); end

%% MAKE OUTPUT STRUCTURE
S = [] ;
for ii = 1:Nsens
	if ii ~= Iz; S = setfield(S,colname{ii},A{ii}(I)) ; end
end
S.dnum = dnum(I) ;

%% BIN AVERAGE
xb   = [0:bsz:max(z(I))+0.5*bsz] ;
flds = fieldnames(S) ;
for ii = 1:Nsens
	yb = bin_data(z(I),getfield(S,flds{ii}),xb) ;	
	S  = setfield(S,flds{ii},yb) ;	
end
S.z  = xb ;

%% EXTRACT INFO FROM FILE NAME
if exist('delim') 
	str = split(fname{kk}(1:end-4),delim); 
	for ii = 1:numel(namef)
		S = setfield(S,namef{ii},str{field{ii}}) ;
	end
end

%% ATTACH GPS COORDINATES
if exist('lon')
	S.lon = interp1(tim,lon,S.dnum) ;
	S.lat = interp1(tim,lat,S.dnum) ; end

%% USE GSW IF POSSIBLE
if gsw
	S.asal	 = gsw_SA_from_SP(S.sal,0,-66.4,49);
	S.ctemp	 = gsw_CT_from_t(S.asal,S.temp,0);
	S.pden	 = gsw_sigma0(S.asal,S.ctemp); end

%% SHOW KEPT/DISCARDED DATA 
if startsWith(showz,'y')
	plot(S.dnum,S.z,'.');
	datetick('x','HHMM');
	ylabel('z (m)');
	title(fname{kk},'interpreter','none') ; end

%% OUTPUT
if nargout<1 
	save([fpath{kk} '/' fname{kk}(1:end-4) '.mat'],'S') ;
elseif nargout==1 & numel(fname)==1
	varargout{1} = S;
end
 
end

end
