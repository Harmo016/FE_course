function [ bsPrice ] = blackScholes( S,K,t,r_f,sigma,type )


d1 = 1/(sigma*sqrt(t))*(log(S/K)+(r_f+(sigma^2)/2)*t);

d2 = d1-sigma*sqrt(t);


switch type
    
    case 'Put'
        bsPrice = normcdf(-d2) * K * exp(-r_f*t) - normcdf(-d1) * S;
        
    case 'Call'
        bsPrice = normcdf(d1) * S - normcdf(d2) * K * exp(-r_f*t);
end

end

