function [C] = optionPricing(S0, K, r, sigma, T, optionType)
% Denna funktion beräknar det analytiska uttrycket enligt givna uttryck.
% optionType = 1 => Cash-or-nothoing-option
% optionType = 2 => Asset-or-nothoing-option

d_1 = (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
d_2 = d_1 - sigma*sqrt(T);

if optionType == 1
    
    C = exp(-r*T)*normcdf(d_2);
    
elseif optionType == 2
    C = S0*normcdf(d_1);

end 


end

