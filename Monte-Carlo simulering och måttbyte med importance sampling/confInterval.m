function [CI, length] = confInterval(lvl,f_hat,sigma_hat,N)
%CONFINTERVAL Summary of this function goes here
%   Detailed explanation goes here
low = f_hat - norminv(1-lvl/2) * sigma_hat / sqrt(N);
upper = f_hat + norminv(1-lvl/2) * sigma_hat / sqrt(N);
length = abs(upper-low);
CI = [low ; upper];
end

