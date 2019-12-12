clear
clc

%% Uppgift 1
S0=1610;
K=1610;
r=0;
sigma=0.163;
mu=0.0166;
T=days252bus('9/10/19','11/15/19')/252;
N_x=5;
N_t=5;

%% Uppgift 1
[S, F] = finiteDifferenceMethods(S0,K,r,sigma,mu,T,N_x,N_t,0,0,1);

[minVal,index] = min(abs(S-S0));

% for i = 1:N_S
% C(i)=blsprice(S(i),K,r,T,sigma);
% end

%Plot med ökande N_S (rummet)

plotS = zeros(20,1);

for i = 1:20
    
    [S, f] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_t+i, N_t,0,0,1);
    
    [minVal,index] = min(abs(S-S0));
    
    plotS(i) = f(index);
    
end

plotT = zeros(20,1);

for i = 1:20
    
    [S, f] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x, N_t+i,0,0,1);
    
    [minVal,index] = min(abs(S-S0));
    
    plotT(i) = f(index);
    
end



%plot(plotS)
%plot(plotT)

%error = C(:)-f(:,1);
deltaLinear = (f(index,1)-f(index+1,1))/(S(index)-S(index+1));
deltaBs = blsdelta(S0,K,r,T,sigma);


%% Uppgift 2

[X_theta_0, grid_theta_0] = finiteDifferenceMethods(S0,K,r,sigma,mu,T,30,30,0,1,2);
[X_theta_05, grid_theta_05] = finiteDifferenceMethods(S0,K,r,sigma,mu,T,30,30,0.5,1,2);
[X_theta_1, grid_theta_1] = finiteDifferenceMethods(S0,K,r,sigma,mu,T,30,30,1,1,2);


[minVal,index_theta_0] = min(abs(exp(X_theta_0)-S0));
[minVal,index_theta_05] = min(abs(exp(X_theta_05)-S0));
[minVal,index_theta_1] = min(abs(exp(X_theta_1)-S0));

%Optionens värde enligt theta = 0
grid_theta_0(index_theta_0,1);
%Optionens värde enligt theta = 0.5
grid_theta_05(index_theta_05,1);
%Optionens värde enligt theta = 1
grid_theta_1(index_theta_1,1);



%Testar stabilitet
%Testar gridstorlek 30x30, 50x50, 100x100, 500x500 för theta = 0, 0.5 och 1
% %Undersöker implicit:
% [X_theta_0_30_30, grid_theta_0_30_30] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,30,30,0,1); %stabil
% [X_theta_0_50_50, grid_theta_0_50_50] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,50,50,0,1); %stabil
% [X_theta_0_100_100, grid_theta_0_100_100] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,100,100,0,1);  %stabil
% [X_theta_0_500_500, grid_theta_0_500_500] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,500,500,0,1); %stabil

%=> Stabil för alla deltaT och deltaX 

%Undersöker Crank-Nicholson
% [X_theta_05_30_30, grid_theta_05_30_30] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,30,30,0.5,1); %stabil
% [X_theta_05_50_50, grid_theta_05_50_50] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,50,50,0.5,1); %stabil
% [X_theta_05_100_100, grid_theta_05_100_100] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,100,100,0.5,1);  %stabil
% [X_theta_05_500_500, grid_theta_05_500_500] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,500,500,0.5,1); %stabil

%=> Stabil för alla deltaT och deltaX 

%Undersöker explicit
% [X_theta_1_30_30, grid_theta_1_30_30] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,30,30,1,1,2); %stabil
% [X_theta_1_35_35, grid_theta_1_35_35] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,35,35,1,1,2); %stabil
% [X_theta_1_40_40, grid_theta_1_40_40] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,40,40,1,1,2); %stabil
% [X_theta_1_45_45, grid_theta_1_45_45] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,45,45,1,1,2); %stabil
% [X_theta_1_50_50, grid_theta_1_50_50] = finiteDifferenceMethods(S0,K,r,sigma,mu,10,50,50,1,1,2); %instabil



%=> deltaT blir ostabil i intervallet [0.2041,0.2273]
%=> deltaX blir ostabil i intervallet [0.0692,0.0771]


% 
%     [S, f] = grid(S0, K, r, sigma, mu, T, N_S, N_t+i);
%     [minVal,index] = min(abs(S-S0));
%     plotT(i) = f(index);


% Undersök tidskonvergens för explicit
absolute_error = zeros(20,1);
for i = 1:20
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x, N_t+i,1,1,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    absolute_error(i) = abs(f_test_1(index)-blsprice(S0,K,r,T,sigma)); 
    
end

% plot(absolute_error)
% xlabel('#deltaT')
% ylabel('Absolute error')
% title('Tidskonvergens explicit FDM')

%%

% Undersök tidskonvergens för Crank-Nicholson
absolute_error = zeros(20,1);
for i = 1:20
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x, N_t+i,0.5,1,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    absolute_error(i) = abs(f_test_1(index)-blsprice(S0,K,r,T,sigma)); 
    
end

% plot(absolute_error)
% xlabel('#deltaT')
% ylabel('Absolute error')
% title('Tidskonvergens C-N')

%%
% Undersök tidskonvergens för implicit
absolute_error = zeros(20,1);
for i = 1:20
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x, N_t+i,0,1,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    absolute_error(i) = abs(f_test_1(index)-blsprice(S0,K,r,T,sigma)); 
    
end

% plot(absolute_error)
% xlabel('#deltaX')
% ylabel('Absolute error')
% title('Tidskonvergens implicit')
 
%%

% Undersök rumskonvergens för explicit
absolute_error = zeros(20,1);
for i = 1:20
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x+i, N_t,1,1,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    absolute_error(i) = abs(f_test_1(index)-blsprice(S0,K,r,T,sigma)); 
    
end

% plot(absolute_error)
% xlabel('#deltaX')
% ylabel('Absolute error')
% title('Rumskonvergens explicit')
%%

% Undersök rumskonvergens för Crank-Nicholson
absolute_error = zeros(20,1);
for i = 1:20
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x+i, N_t,0.5,1,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    absolute_error(i) = abs(f_test_1(index)-blsprice(S0,K,r,T,sigma)); 
    
end

% plot(absolute_error)
% xlabel('#deltaX')
% ylabel('Absolute error')
% title('Rumskonvergens C-N')

%%
% Undersök rumskonvergens för implicit
absolute_error = zeros(20,1);
for i = 1:20
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x+i, N_t,0,1,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    absolute_error(i) = abs(f_test_1(index)-blsprice(S0,K,r,T,sigma)); 
    
end

% plot(absolute_error)
% xlabel('#deltaX')
% ylabel('Absolute error')
% title('rumskonvergens implicit')



%% uppgift 3

%Pris för cash-or-nothing-option
[con_S0, con_grid]=finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x, N_t,0,2,2);
[minVal,index] = min(abs(exp(con_S0)-S0));
C_con =con_grid(index);

%Pris för asset-or-nothing-option
[aon_S0, aon_grid]=finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x, N_t,0,3,2);
[minVal,index] = min(abs(exp(aon_S0)-S0));
C_aon =aon_grid(index);


%% Undersöker relativa felet för alla 3 optioner med FDM-metoderna implicit samt Crank-Nicholson
% Undersök rumskonvergens för explicit


%vanlig option, implicit
relative_error_call_implicit = zeros(3,1);
for i = 1:3
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x+i*10, N_t+i*10,0,1,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    relative_error_call_implicit(i) = f_test_1(index)/blsprice(S0,K,r,T,sigma)-1; 
    
end


%vanlig option, Crank
relative_error_call_Crank = zeros(3,1);
for i = 1:3
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x+i*10, N_t+i*10,0.5,1,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    relative_error_call_Crank(i) = f_test_1(index)/blsprice(S0,K,r,T,sigma)-1; 
    
end



%cash-or-nothing, implicit
relative_error_con_implicit = zeros(3,1);
for i = 1:3
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x+i*10, N_t+i*10,0,1,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    relative_error_con_implicit(i) = f_test_1(index)/optionPricing(S0,K,r,T,sigma,1)-1; 
    
end


%cash-or-nothing, Crank
relative_error_con_Crank = zeros(3,1);
for i = 1:3
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x+i*10, N_t+i*10,0.5,2,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    relative_error_con_Crank(i) = f_test_1(index)/optionPricing(S0,K,r,T,sigma,1)-1; 
    
end


%asset-or-nothing, implicit
relative_error_aon_implicit = zeros(3,1);
for i = 1:3
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x+i*10, N_t+i*10,0,3,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    relative_error_aon_implicit(i) = f_test_1(index)/optionPricing(S0,K,r,T,sigma,2)-1; 
    
end
%plot(relative_error_aon_implicit);
%hold on

%asset-or-nothing, Crank
relative_error_aon_Crank = zeros(3,1);
for i = 1:3
    
    [S_test_1, f_test_1] = finiteDifferenceMethods(S0, K, r, sigma, mu, T, N_x+i*10, N_t+i*10,0.5,3,2);
    [minVal,index] = min(abs(exp(S_test_1)-S0));
    relative_error_aon_Crank(i) = f_test_1(index)/optionPricing(S0,K,r,T,sigma,2)-1; 
    
end


%plot(relative_error_aon_Crank);


