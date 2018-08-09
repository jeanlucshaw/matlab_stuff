function m_circle(lon,lat,r,label)

if isempty(label)
	rng = [0:360 0] ;
else
	rng = [290:360 0:250] ;
	m_text(lon,lat-r,label,'horizontalalignment','center','fontsize',6) ;
end

x     = lon + r*cosd(rng)/cosd(lat) ;
y     = lat + r*sind(rng) ;
m_plot(x,y,'k-') ;

end
