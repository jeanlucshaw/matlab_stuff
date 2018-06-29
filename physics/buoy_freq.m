function N = buoy_freq(z,rho) 

g       = 9.81 ;
rho0    = nanmean(rho) ;

rho     = sort(rho) ;

kmax    = length(z) ;

for k = 1:kmax
    
    if     k == 1
        dz   = z(k+1)    -    z(k);
        drho = rho(k+1)  -    rho(k) ;
    elseif k == kmax
        dz   = z(k)      -    z(k-1);
        drho = rho(k)    -    rho(k-1) ;
    else
        dz   = z(k+1)    -    z(k-1);
        drho = rho(k+1)  -    rho(k-1) ;    
    end

    N(k) = sqrt((g/rho0)*drho/dz) ;

end

% oubedon  N = sqrt((g/rho0)*gradient(rho,z))