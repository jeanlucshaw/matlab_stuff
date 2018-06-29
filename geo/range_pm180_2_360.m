function output = range_pm180_2_360(input)

I = find(input < 0) ;
input(I) = input(I) + 360 ;

output = input ;

end
