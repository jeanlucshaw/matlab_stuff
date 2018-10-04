function [time,val]	= monthly(intime,inval,mode)
% Syntax : 	[time,val]	= monthly(intime,inval,mode)
%
% Monthly statistics. Depending on the the value of 'mode', this routine
% returns the monthly mean, standard deviation, maximum or minimum of vector
% 'inval' corresponding to times 'intime'. Time is expected in matlab format
% as days since january 1, year 0. It may bug for years before 1000.
%
% Mode is expected to be a string array, either 'mean', 'std', 'max', or 'min'.  


	% Sort
	[t,I]	= sort(intime) ;
	v	= inval(I) ;	

	% Find/last first year and first month
	fy	= str2num(datestr(t(1)  ,'yyyy')) ;  
	fm	= str2num(datestr(t(1)  ,'mm')) ;  
	ly	= str2num(datestr(t(end),'yyyy')) ;  
	lm	= str2num(datestr(t(end),'mm')) ;  

	% Init output
	val	= [] ;
	time 	= [] ;	
	for ii 	= fy:ly
		if ii == fy & ii == ly
			FM	= fm ;
			LM	= lm ;
		elseif 	ii == fy
			FM	= fm ;
			LM	= 12 ;
		elseif	ii == ly
			FM	= 1  ;
			LM	= lm ;
		else	
			FM	= 1  ;
			LM	= 12 ;
		end

		for jj	= FM:LM
			as	= datenum([sprintf('%04d',ii) sprintf('%02d',jj)],'yyyymm') ;
			if jj == 12
				af	= datenum([sprintf('%04d',ii+1) '01'                  ],'yyyymm') ;
			else
				af	= datenum([sprintf('%04d',ii  ) sprintf('%02d',jj+1)  ],'yyyymm') ;
			end
			
			I	= find( t < af  &  t > as ) ;
			switch mode
				case 'mean'
					nuv	= mean(v(I),'omitnan') ;
					nut	= mean(t(I),'omitnan') ;
				case 'std'
					nuv	=  std(v(I),'omitnan') ;
					nut	=  std(t(I),'omitnan') ;
				case 'max'
					[nuv,II]	= max(v(I),[],'omitnan') ;
					nut		= t(II) ;
				case 'min'
					[nuv,II]	= min(v(I),[],'omitnan') ;
					nut		= t(II) ;
			end

			val	= app(val ,nuv) ;
			time	= app(time,nut) ;
		end
	end
end
