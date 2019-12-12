function [X, grid] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x, N_t, theta,optionType,valueMethod)
% Denna funktion beräknar numeriska optionsvärden för 3 olika optioner med
% två olika tillvägagångssätt.

% Option type 1 = vanlig köpoption
% Option type 2 = cash-or-nothing-option
% Option type 3 = asset-or-nothing-option
% Value method 1 = Hulls metod
% Value method 2 = Andersen och Brotherton-Ratcliff metod


if valueMethod == 1
    
    
    a=log(S0)+(mu-(sigma^2/2))*T;
    b=sigma^2*T;
    k=norminv(0.9995);
    
    S_lower=exp(a-k*sqrt(b));
    S_upper=exp(a+k*sqrt(b));
    
    delta_t=T/(N_t-1)
    delta_S=(S_upper-S_lower)/(N_x-1)
    
    grid=zeros(N_x,N_t);
    
    X=zeros(N_x,1);
    
    
    
    for i = 0:(N_x-1)
        X(i+1)=S_upper - delta_S*i;
    end
    
    
    for i = 0:(N_t-1)
        grid(1,i+1)= max(S_upper-K*exp(-r*(T-delta_t*i)),0);
    end
    
    for i = 0:(N_t-1)
        grid(N_x,i+1)= max(S_lower-K*exp(-r*(T-delta_t*i)),0);
    end
    
    for i = 1:N_x
        grid(i,N_t)= max(S_upper-(i-1)*delta_S-K,0);
    end
    
    
    for j = (N_t-1):-1:1
        for i = 2:(N_x-1)
            grid(i,j) = (1/(r+(1/delta_t)))*((grid(i,j+1)/delta_t)+r*X(i)*((grid(i,j+1)-grid(i+1,j+1))/(2*delta_S))+0.5*(sigma^2)*(X(i)^2)*((grid(i-1,j+1)-2*grid(i,j+1)+grid(i+1,j+1))/(delta_S^2)));
        end
    end
    
    
    
end


if valueMethod == 2
    
    
    
    X0 = log(S0);
    a=log(S0)+(mu-(sigma^2/2))*T;
    b=sigma^2*T;
    k=norminv(0.9995);
    cash = 1;
    asset = S0;
    
    X_lower=a-k*sqrt(b);
    X_upper=a+k*sqrt(b);
    
    delta_t=T/(N_t-1);
    delta_x=(X_upper-X_lower)/(N_x-1);
    
    H=zeros(N_x,N_t);
    
    X=zeros(N_x,1);



for i = 0:(N_x-1)
    X(i+1)=X_upper - delta_x*i;
end

if optionType == 1
    
    
    for i = 0:(N_t-1)
        H(1,i+1)= max(exp(X_upper)-K*exp(-r*(T-delta_t*i)),0);
    end
    
    for i = 0:(N_t-1)
        H(N_x,i+1)= max(exp(X_lower)-K*exp(-r*(T-delta_t*i)),0);
    end
    
    for i = 1:N_x
        H(i,N_t)= max(exp(X_upper-(i-1)*(delta_x))-K,0);
    end
    
    
elseif optionType == 2
    
    for i = 0:(N_t-1)
        H(1,i+1)= max(cash*exp(-r*(T-delta_t*i)),0);
    end
    
    for i = 0:(N_t-1)
        if exp(X_lower)-K*exp(-r*(T-delta_t*i))<0
            H(N_x,i+1)= 0;
        else
            H(N_x,i+1)= cash*exp(-r*(T-delta_t*i));
        end
    end
    
    for i = 1:N_x
        if exp(X_upper-(i-1)*delta_x)>K
            H(i,N_t)=cash;
        end
    end
    
elseif optionType == 3
    
    for i = 0:(N_t-1)
        H(1,i+1)= max(asset*exp(-r*(T-delta_t*i)),0);
    end
    
    for i = 0:(N_t-1)
        if exp(X_lower)-K*exp(-r*(T-delta_t*i))<0
            H(N_x,i+1)= 0;
        else
            H(N_x,i+1)= asset*exp(-r*(T-delta_t*i));
        end
    end
    
    
    for i = 1:N_x
        if exp(X_upper-(i-1)*delta_x)>K
            H(i,N_t)=asset;
        end
    end
    
end

c = -delta_t/(delta_x^2)*sigma^2;
u = 0.5*delta_t/(delta_x^2)*(sigma^2+delta_x*(r-0.5*sigma^2));
l = 0.5*delta_t/(delta_x^2)*(sigma^2-delta_x*(r-0.5*sigma^2));

M= zeros(N_x-2,N_x-2);
M(1,1) = c;
M(1,2) = u;
M(N_x-2,N_x-3) = l;
M(N_x-2,N_x-2) = c;



for i = 2:(N_x-3)
    M(i, i-1) = l;
    
end

for i = 2:(N_x-3)
    M(i, i) = c;
end

for i = 2:(N_x-3)
    M(i, i+1) = u;
end

I=eye(N_x-2,N_x-2);

B = zeros(N_x-2,N_t-1);
for j= (N_t-1):-1:1
    B(1,j) = l*((1-theta)*H(1,j)+theta*H(1,j+1));
    B(N_x-2,j) = u*((1-theta)*H(N_x,j)+theta*H(N_x,j+1));
end

B(1,N_t-1) = l*((1-theta)*H(1,N_t-1)+theta*H(1,N_t));
B(N_x-2,N_t-1) = u*((1-theta)*H(N_x,N_t-1)+theta*H(N_x,N_t));


stora_H_rand = H(2:(N_x-1),:);
stora_H = zeros(N_x-2,N_t-1);


for i = (N_t-1):-1:1
    stora_H(:,i) =  inv((1+r*delta_t).*I-(1-theta).*M)*((theta.*M+I)*stora_H_rand(:,i+1)+B(:,i));
    stora_H_rand(:,i) = stora_H(:,i);
end

grid = zeros(N_x,N_t);

grid(2:N_x-1,:) = stora_H_rand;
grid(1,:) = H(1,:);
grid(N_x,:) = H(N_x,:);


x_check=exp(X);

end 



