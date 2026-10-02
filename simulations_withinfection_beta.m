close all;
clear all;

%% --- Base parameters (north) ---
Lambda = [300; 50; 659.8; 27; 36.1; 0.32];

mu = zeros(length(Lambda),1);
b_hat_zL = [73.4; 44.5; 87.8; 165.4; 28.9; 1963.3];
b_hat_zN = [3.7; 34.2; 0.35; 277.3; 0.36; 750].*[1; 1; 3.59; 1; 1.07; 1.25];

base_c = 0.1815;
b_hat_L = b_hat_zL./[base_c; base_c; base_c; base_c; 1-base_c; base_c];
b_hat_N = b_hat_zN./[base_c; base_c; base_c; base_c; 1-base_c; base_c];
bA = 1037.5*1.08;

betaL = [0.91; 0.57; 0.56; 0.06; 0.024];
betaN = [0.91; 0.57; 0.56; 0.06; 0.024];
alpha = [1;1;1;1;1];

k = [2; 2; 2; 2; 1];
r = 2.52e6; 
K = 5;
sL = 0.43; sN = 0.55; sA = 0.84;
n_hosts = length(Lambda);

% %% --- Base parameters (south) ---
% Lambda = [273.5; 50; 174.9; 27; 1278.1; 0.21];
% 
% mu = zeros(length(Lambda),1);
% b_hat_zL = [3.1; 4.5; 34.5; 16.54; 25.2; 196.3];
% b_hat_zN = [0.25; 3.42; 0.035; 27.73; 2.5; 75.0].*[1.24; 1; 0.90; 1; 1.19; 1.00];
% 
% base_c = 0.0204;
% b_hat_L = b_hat_zL./[base_c; base_c; base_c; base_c; 1-base_c; base_c];
% b_hat_N = b_hat_zN./[base_c; base_c; base_c; base_c; 1-base_c; base_c];
% bA = 1037.5*1.00;
% 
% betaL = [0.91; 0.57; 0.56; 0.06; 0.24];
% betaN = [0.91; 0.57; 0.56; 0.06; 0.24];
% alpha = [1;1;1;1;1];
% 
% sL = 0.43; sN = 0.55; sA = 0.84;
% k = [2; 2; 2; 2; 1];
% n_hosts = length(Lambda);
% r = 2.52e6; K=5;

% Time simulation
T = 1000;
beta_values = linspace(0,5,100);
n_c = length(beta_values);

%% Population Vectors & Functions

L = zeros(1,T+1); N = zeros(1,T+1); A = zeros(1,T+1);
LS = zeros(1,T+1); NS = zeros(1,T+1); AS = zeros(1,T+1);
LI = zeros(1,T+1); NI = zeros(1,T+1); AI = zeros(1,T+1);
H = zeros(n_hosts,T+1); HS = zeros(n_hosts,T+1); HI = zeros(n_hosts,T+1);

pL = zeros(n_hosts); pN = zeros(n_hosts);
SL = zeros(1,T); SN = zeros(1,T); SA = zeros(1,T);
infect_exponential = zeros(n_hosts-1, 1);

% Preallocate storage for final results across c
L_final = zeros(1,n_c); N_final = zeros(1,n_c); A_final = zeros(1,n_c);
LI_final = zeros(1,n_c); NI_final = zeros(1,n_c); AI_final = zeros(1,n_c);
LS_final = zeros(1,n_c); NS_final = zeros(1,n_c); AS_final = zeros(1,n_c);
H_final = zeros(n_hosts,n_c);

% %% Initital ConditionsA0
%North
L(1) = 115000; N(1) = 11500; A(1) = 3450;
L(1) = L(1); N(1) = N(1); A(1) = A(1);
LIP(1) = 0; NIP(1) = 0.15; AIP(1) = 0.15;
LS(1) = L(1); LI(1) = 0; L(1) = L(1);
NS(1) = N(1)-NIP(1)*N(1); NI(1) = NIP(1)*N(1); N(1) = N(1);
AS(1) = A(1)-AIP(1)*A(1); AI(1) = A(1)*AIP(1); A(1) = A(1);

H(:,1) = [40; 30; 59.2; 6; 27; 0.46];
HS(:,1) = H(:,1) - [0.15; 0.15; 0.15; 0.15; 0.05; 0].*H(:,1);
HI(:,1) = [0.15; 0.15; 0.15; 0.15; 0.05; 0].*H(:,1);

% %South
% L(1) = 115000; N(1) = 11500; A(1) = 3450;
% L(1) = L(1); N(1) = N(1); A(1) = A(1);
% LIP(1) = 0; NIP(1) = 0.015; AIP(1) = 0.015;
% LS(1) = L(1); LI(1) = 0; L(1) = L(1);
% NS(1) = N(1)-NIP(1)*N(1); NI(1) = NIP(1)*N(1); N(1) = N(1);
% AS(1) = A(1)-AIP(1)*A(1); AI(1) = A(1)*AIP(1); A(1) = A(1);
% 
% H(:,1) = [27.885; 30; 27.742; 6; 186.994; 0.40];
% HS(:,1) = H(:,1) - [0.015; 0.015; 0.015; 0.015; 0.005; 0].*H(:,1);
% HI(:,1) = [0.015; 0.015; 0.015; 0.015; 0.005; 0].*H(:,1);

%% Compute mu and zs

for i = 1:n_hosts-1
    mu(i) = k(i)*log(Lambda(i)/(k(i)*H(i,1))+1);
end
mu(n_hosts) = 2*log((Lambda(n_hosts)+sqrt(Lambda(n_hosts)^2+4*H(n_hosts,1)^2))/(2*H(n_hosts,1)));


SLsum0 = 0; SNsum0 = 0; 
for i = 1:n_hosts-1
    SLsum0 = SLsum0+b_hat_zL(i)*(H(i,1)+Lambda(i)/k(i));
    SNsum0 = SNsum0+b_hat_zN(i)*H(i,1);
end
SLsum0 = SLsum0+b_hat_zL(n_hosts)*(H(n_hosts,1)*exp(-mu(n_hosts)/2)+Lambda(n_hosts));
SNsum0 = SNsum0+b_hat_zN(n_hosts)*H(n_hosts,1)*exp(-mu(n_hosts)/2);
SAsum0 = bA*H(n_hosts,1);

zL = -L(1)/sL*log(1-0.1/sL)/SLsum0;
zN = -N(1)/sN*log(1-0.1/sN)/SNsum0; SN(1) = 1-exp(-zN/(N(1)/sN)*SNsum0);
zA = -sN*SN(1)*(N(1)/sN)*log(1-0.3/sA)/SAsum0;

bL = zeros(n_hosts,1); bN = zeros(n_hosts,1);
idx = [1:(n_hosts-2) n_hosts];
for i = idx
    bL(i) = b_hat_L(i) * base_c;
    bN(i) = b_hat_N(i) * base_c;
end
bL(n_hosts-1) = b_hat_L(n_hosts-1) * (1-base_c);
bN(n_hosts-1) = b_hat_N(n_hosts-1) * (1-base_c);

for j=1:n_c
    beta_new = beta_values(j)*ones(5,1);
    L(1) = 115000; N(1) = 11500; A(1) = 3450;
    LS(1) = L(1); LI(1) = 0;
    NS(1) = N(1)-NIP(1)*N(1); NI(1) = NIP(1)*N(1);
    AS(1) = A(1)-AIP(1)*A(1); AI(1) = A(1)*AIP(1);
    %% Time loop
    for t = 1:T
        % pL and pN
        pL_denom = 0; pN_denom = 0;
        
        for i = 1:n_hosts-1          
            pL_denom = pL_denom+bL(i)*(H(i,t)+Lambda(i)/k(i))*exp(-mu(i)/2);
            pN_denom = pN_denom+bN(i)*H(i,t);
        end
        pL_denom = pL_denom+bL(n_hosts)*H(n_hosts,t);
        pN_denom = pN_denom+bN(n_hosts)*H(n_hosts,t)*exp(-mu(n_hosts)/2);
    
        for i = 1:n_hosts-1
            pL(i) = bL(i)*(H(i,t)+Lambda(i)/k(i))*exp(-mu(i)/2)/pL_denom;
            pN(i) = bN(i)*H(i,t)/pN_denom;
        end
        pL(n_hosts) = bL(n_hosts)*H(n_hosts,t)/pL_denom;
        pN(n_hosts) = bN(n_hosts)*H(n_hosts,t)*exp(-mu(n_hosts)/2)/pN_denom;

        % SLsum, SNsum, SAsum
        SLsum = 0; SNsum = 0; 
        for i = 1:n_hosts-1
            SLsum = SLsum+bL(i)*(H(i,t)+Lambda(i)/k(i));
            SNsum = SNsum+bN(i)*H(i,t);
        end
        SLsum = SLsum+bL(n_hosts)*(H(n_hosts,t)*exp(-mu(n_hosts)/2)+Lambda(n_hosts));
        SNsum = SNsum+bN(n_hosts)*H(n_hosts,t)*exp(-mu(n_hosts)/2);
        SAsum = bA*H(n_hosts,t);

        % birth term
        numerator = r*0.5*A(t); denominator = 0.5*A(t)/H(n_hosts,t)+K;
        birth_term = numerator/denominator;

        % Compute SL, SN, SA
        SL(t) = 1-exp(-zL/birth_term*SLsum);
        SN(t) = 1-exp(-zN/L(t)*SNsum);
        SA(t) = 1-exp(-zA/(sN*SN(t)*L(t))*SAsum);
    
        % infection dynamics
        infectionLsum = 0; infectionNsum = 0; infectionLsum2 = 0;
        for i = 1:n_hosts-1
            infect_exponential(i) = exp(-pN(i)*alpha(i)*sN*SN(t)*LI(t)/H(i,t));
            infectionLsum = infectionLsum+pL(i)*beta_new(i)*(HI(i,t))/(H(i,t)+Lambda(i)/k(i));
            infectionNsum = infectionNsum+pN(i)*beta_new(i)*HI(i,t)*exp(-mu(i))/H(i,t);
        end
    
        % Tick Equations
        L(t+1) = birth_term*sL*SL(t);
        N(t+1) = sN*SN(t)*L(t);
        A(t+1) = sA*SA(t)*sN*SN(t)*L(t);
        LS(t+1) = birth_term*sL*SL(t)*(1-infectionLsum);
        LI(t+1) = birth_term*sL*SL(t)*infectionLsum;
        NS(t+1) = sN*SN(t)*LS(t)*(1-infectionNsum);
        NI(t+1) = sN*SN(t)*(LI(t)+LS(t)*infectionNsum);
        AS(t+1) = sA*SA(t)*NS(t+1);
        AI(t+1) = sA*SA(t)*NI(t+1);
    
       % Host Equations
        for i = 1:n_hosts-2
            HS(i,t+1) = (HS(i,t)+Lambda(i)/2*(1+exp(mu(i)/2)))*exp(-mu(i))*infect_exponential(i);
            HI(i,t+1) = (HI(i,t)+(HS(i,t)+Lambda(i)/2*(1+exp(mu(i))))*(1-infect_exponential(i)))*exp(-mu(i));
            H(i,t+1) = (H(i,t)+Lambda(i)/2*(1+exp(mu(i)/2)))*exp(-mu(i));
        end
        HS(n_hosts-1,t+1) = (HS(n_hosts-1,t)+Lambda(n_hosts-1))*exp(-mu(n_hosts-1))*infect_exponential(n_hosts-1);
        HI(n_hosts-1,t+1) = (HI(n_hosts-1,t)+(HS(n_hosts-1,t)+Lambda(n_hosts-1))*(1-infect_exponential(n_hosts-1)))*exp(-mu(n_hosts-1));
        H(n_hosts-1,t+1) = (H(n_hosts-1,t)+Lambda(n_hosts-1))*exp(-mu(n_hosts-1));
    
        H(n_hosts,t+1) = (H(n_hosts,t)*exp(-mu(n_hosts)/2)+Lambda(n_hosts))*exp(-mu(n_hosts)/2);
        HS(n_hosts,t+1) = (H(n_hosts,t)*exp(-mu(n_hosts)/2)+Lambda(n_hosts))*exp(-mu(n_hosts)/2);
        HI(n_hosts,t+1) = 0;
    end

    for i = 1:n_hosts-2
            HS(i,t+1) = (HS(i,t)*infect_exponential(i)+Lambda(i)/2*(1+exp(mu(i)/2)))*exp(-mu(i));
            HI(i,t+1) = (H(i,t)-HS(i,t)*infect_exponential(i))*exp(-mu(i));
            H(i,t+1) = (H(i,t)+Lambda(i)/2*(1+exp(mu(i)/2)))*exp(-mu(i));
    end
    
    HS0 = mean(HS(:,end-10:end),2);
    HI0 = mean(HI(:,end-10:end),2);
    for i=1:n_hosts
        HI_final(i,j) = mean(HI(i,end-10:end));
    end
    for i=1:n_hosts
        HS_final(i,j) = mean(HS(i,end-10:end));
    end
    for i=1:n_hosts
        H_final(i,j) = mean(H(i,end-10:end));
    end

end
figure;
linestyles = {'-', '--', ':', '-.'};

for i=1:n_hosts-1
    styleIdx = mod(i-1, length(linestyles)) + 1;
    plot(beta_values,HI_final(i,:), 'LineStyle',linestyles{styleIdx}, 'LineWidth', 2);
    hold on
end
xlabel('$\beta$', 'FontSize', 14, 'Interpreter', 'latex'); ylabel('$H_{i_I}$ (host/ha)', 'FontSize', 14, 'Interpreter', 'latex');
legend('Mice','Chipmunks','Shrews','Squirrels','Lizards')
grid on;
figure;
for i=1:n_hosts
    plot(beta_values,HS_final(i,:), 'LineWidth', 2);
    hold on
end
xlabel('beta'); ylabel('Susceptible Host');
legend('Mice','Chipmunks','Shrews','Squirrels','Lizards','Deer')
grid on;
figure;
for i=1:n_hosts
    plot(beta_values,H_final(i,:), 'LineWidth', 2);
    hold on
end
xlabel('beta'); ylabel('Total Host');
legend('Mice','Chipmunks','Shrews','Squirrels','Lizards','Deer')
grid on;