function output = range_360_2_pm180(input)

I = find(input > 180) ;
input(I) = input(I) - 360 ;

output = input ;

end
