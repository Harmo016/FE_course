function [S , xi] = importanceSampling(S_0,sigma,rf,M,T,phi)
%BSM Summary of this function goes here
%   Detailed explanation goes here
deltaT=T/M;
S = zeros(M+1,1);
S(1) = S_0;

% sätta värde på drift justering
% motivering?



    rnd = normrnd(0,1);
    S(2) = S(1)*exp(((rf-(sigma^2)/2) + phi * sigma)*deltaT+rnd*sigma*sqrt(deltaT));

    g_x = normpdf(rnd+phi * deltaT,phi * deltaT,1);
    f_x = normpdf(rnd+phi * deltaT);
    xi = f_x/g_x;

end

