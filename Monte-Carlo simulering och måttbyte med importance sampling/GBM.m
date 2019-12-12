function [S] = GBM(S_0,sigma,rf,M,T)
%BSM Summary of this function goes here
%   Detailed explanation goes here
deltaT=T/M;
S = zeros(M+1,1);
S(1) = S_0;
for i=2:M+1
    S(i) = S(i-1)*exp((rf-(sigma^2)/2)*deltaT+normrnd(0,1)*sigma*sqrt(deltaT));
end
end

