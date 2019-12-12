% TPPE53 - Finansiella värderingsmetodiker ----- Inlämning 2
% Erik Karvonen - erika797
% Harald Graf Morin - harm016
%%------------------------------------------------------------------------
clc
clear all
close all
%% Uppgift 1

N = 1E4; %Skala med k=4 för att halvera längden på intervallet
mu = 0;
sigma = 0.2;

%% Implementera Monte Carlo-simulering utan variansreducering
% Simulera X

X = normrnd(mu,sigma,N,1);
f = exp(X);

% Skatta f
f_hat = mean(exp(X));

% Skatta std

sigma_hat = sqrt((1/(N-1))*sum((f-f_hat).^2));

% Konfidensintervall

conf_lvl = 0.05;

[CI_f_MC, length_MC] = confInterval(conf_lvl,f_hat,sigma_hat,N);


%% Antitetisk 
% Skapa likformiga par
% Nytt f av paren (medel)
% kötta på som vanligt

U_vector = rand(N,1); % rand genererar likformigt fördelade slumptal
U_vector_mirr = 1-U_vector;

X_a = norminv(U_vector, mu, sigma);
X_a_mirr = norminv(U_vector_mirr, mu, sigma);

f_a_u = exp(X_a);
f_a_tilde = exp(X_a_mirr);
f_a = (f_a_u+f_a_tilde)/2;

%Skattning medelvärdet
f_hat_a = mean(f_a);
%Skatta stdasiatPO

sigma_hat_a = sqrt((1/(N-1))*sum((f_a-f_hat_a).^2));
% Konfidensintervall
[CI_f_a,length_a] = confInterval(conf_lvl,f_hat_a,sigma_hat_a,N);

%% Latin hyperkub
%Generera slumptal
X_l = lhsnorm(mu,sigma,N);

f_l = exp(X_l);
f_hat_l = mean(f_l);

%skatta std
sigma_hat_l = sqrt((1/(N-1))*sum((f_l-f_hat_l).^2));

[CI_f_l,length_l] = confInterval(conf_lvl,f_hat_l,sigma_hat_l,N);

%% Test konf
M = 1000; %nr of tests
X_est = zeros(N,M);
X_est_var = zeros(N,M);
f_est_var = zeros(N,M);
for i=1:M
    X_est(:,i) = normrnd(mu,sigma,N,1);
    
    U_vector_tmp = rand(N,1); % rand genererar likformigt fördelade slumptal
    U_vector_mirr_tmp = 1-U_vector_tmp;
    
    X_tmp = norminv(U_vector_tmp, mu, sigma);
    X_tmp_tilde = norminv(U_vector_mirr, mu, sigma);
   
    f_est_var(:,i) = (exp(X_tmp)+exp(X_tmp_tilde))/2;
end

f_est = exp(X_est);

f_hat_est = mean(f_est,1);
f_hat_est_var = mean(f_est_var,1);


% Skatta std

sigma_hat_est = sqrt((1/(N-1))*sum((f_est-f_hat_est).^2));
sigma_hat_est_var = sqrt((1/(N-1))*sum((f_est_var-f_hat_est_var).^2));

% Konfidensintervall, antal test kolumner, rad 1 = lower, rad 2 = upper
CI_est = confInterval(conf_lvl,f_hat_est,sigma_hat_est,N);
CI_est_var = confInterval(conf_lvl,f_hat_est_var,sigma_hat_est_var,N);


figure
subplot(2,1,1)
plot(CI_est(1,:))
hold on
plot(CI_est(2,:))
%hold on
%plot(f_est(1,:))

subplot(2,1,2)
plot(CI_est_var(1,:))
hold on
plot(CI_est_var(2,:))
% hold on
% plot(f_est_var(1,:))


% för varje test, kolla vilka som ligger inom intervallet ja => +1. kolla
% andelen sedan med att procenten av ja/alla
true_f=exp((sigma^2)/2);
j_est=zeros(1,M);
for i = 1:M %test
        if true_f >= CI_est(1,i)
            if true_f <= CI_est(2,i)
                j_est(1,i)=j_est(1,i)+1;
            end
        end
end

hits_est = sum(j_est)/M;
j_est_var=zeros(1,M);
for i = 1:M %test
        if true_f >= CI_est_var(1,i)
            if true_f <= CI_est_var(2,i)
                j_est_var(i)=j_est_var(i)+1;
            end
        end
end

hits_est_var = sum(j_est_var)/M;


%% Uppgift 2 - värdera asiatisk option
%Data optionen
T= 1;
K = 100;
vol_y = 0.2;
vol = (1/sqrt(12))* vol_y;
r_f = 0.01;
M=12;

S_0 = 100;
% Monte carlo simulera S_t t=0,1,...,12 
% Nst scenarier - 10 000
N = 1E5;

%loopa Ngr

asiansPOs = zeros(N,1);
Z = zeros(N,1);
for i =1:N

    S = GBM(S_0,vol_y,r_f,M,T);
    asiatPO = max(mean(S)-K,0);
    asiansPOs(i,1)= asiatPO;
    Z(i) = sqrt(S(7)*S(13)); %Kontrollvariat-metoden

end

asianhat = mean(asiansPOs);

sigma_asianhat = sqrt((1/(N-1))*sum((asiansPOs-asianhat).^2));

conf_lvl = 0.05;

[CI_f_asian] = confInterval(conf_lvl,asianhat,sigma_asianhat,N);

%proxy = asianbyhhm(r_f,stockspec(vol,100),'call',K,0,1);


%%  Kontrolvariat

%Analytisk väntevärde Z
my_A = log(S_0) + 1/2 * (r_f- vol_y^2 / 2)*(0.5 + 1);

vol_A = sqrt(vol_y^2*(0.5 + (1/4)*(1-0.5)));

EZ = exp(my_A+vol_A^2 / 2);


Z_hat = mean(Z);

c = -cov(asiansPOs,Z)/var(Z);
% c = -cov(asiansPOs,Z)/var(Z);

asiansPOs_CV = zeros(N,1);
for i =1:N
    asiansPOs_CV(i)= asiansPOs(i) + c(2,1) * (Z(i)-EZ);
end

asianshat_CV = mean (asiansPOs_CV);
sigma_asianhat_CV = sqrt((1/(N-1))*sum((asiansPOs_CV-asianshat_CV).^2));

conf_lvl = 0.05;

[CI_f_asian_CV] = confInterval(conf_lvl,asianshat_CV,sigma_asianhat_CV,N);

%% Uppgift 3 - CoN - Prisestimat utan varians och med importance sampling, 3 sim per metod studera variationen

S_0 = 100;
K = 2 * S_0;
T = 1;
sigma = 0.15;
r_f = 0.01; 
C = 1E8;
N = 1E4;

d2 = 1/(sigma*sqrt(T))*(log(S_0/K)+(r_f-(sigma^2)/2)*T);

CoN_analytic = C * normcdf(d2);
nr_tests = 3;
CoN_hat = zeros(nr_tests,1);
CoN_I_hat = zeros(nr_tests,1);
for j=1:nr_tests
CoN_POs = zeros(N,1);
for i =1:N
    S = GBM(S_0,sigma,r_f,1,T);
    if (S(end))-K>0
       CoN_PO = C;
    else
        CoN_PO = 0;
    end
    CoN_POs(i,1)= CoN_PO;
end

CoN_hat(j) = mean(CoN_POs);

CoN_POs_I = zeros(N,1);

%xi = exp(1+10^2/2);

mu_y = 1/sigma;
for i =1:N
    [S, xi] = importanceSampling(S_0,sigma,r_f,1,T,mu_y);
    if ((S(end))-K) >0
       CoN_PO_I = C * xi;
    else
        CoN_PO_I = 0;
    end
    CoN_POs_I(i,1)= CoN_PO_I;
end

CoN_I_hat(j) = mean(CoN_POs_I);

end

relError_iS = abs(CoN_I_hat - CoN_analytic)./CoN_I_hat;

