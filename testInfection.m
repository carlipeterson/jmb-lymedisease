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

%% Initital Conditions (north)
L0 = 115000;
N0 = 11500;
A0 = 3450;
H(:,1) = [40; 30; 59.2; 6; 27; 0.46];

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
% %% Initital Conditions (south)
% L0 = 115000;
% N0 = 11500;
% A0 = 3450;
% H(:,1) = [27.9; 30; 30; 6; 445.6; 0.40];
%% Parameters to vary
c_values = linspace(0,1,20);

n_c = length(c_values);

%% --- Compute mu and z values ---

for i = 1:n_hosts-1
    mu(i) = k(i)*log(Lambda(i)/(k(i)*H(i,1))+1);
end
mu(n_hosts) = 2*log((Lambda(n_hosts)+sqrt(Lambda(n_hosts)^2+4*H(n_hosts,1)^2)) / (2*H(n_hosts,1)));

SLsum0 = 0; SNsum0 = 0; 
for i = 1:n_hosts-1
    SLsum0 = SLsum0 + b_hat_zL(i)*(H(i,1)+Lambda(i)/k(i));
    SNsum0 = SNsum0 + b_hat_zN(i)*H(i,1);
end
SLsum0 = SLsum0 + b_hat_zL(n_hosts)*(H(n_hosts,1)*exp(-mu(n_hosts)/2)+Lambda(n_hosts));
SNsum0 = SNsum0 + b_hat_zN(n_hosts)*H(n_hosts,1)*exp(-mu(n_hosts)/2);
SAsum0 = bA*H(n_hosts,1);

zL = -L0/sL*log(1-0.1/sL)/SLsum0;
zN = -N0/sN*log(1-0.1/sN)/SNsum0; SN0 = 1-exp(-zN/(N0/sN)*SNsum0);
zA = -sN*SN0*(N0/sN)*log(1-0.3/sA)/SAsum0;

%% Set up for solver and solution storage
unique_solutions1 = []; % will hold all unique soltuions to first gen system
unique_solutions2 = []; % will hold all unique solutions to second gen system

opts = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6); % 

nonmatching_fp2 = cell(n_c,1);      % stores fixed points from system 2 that do not match the system 1 solutions
stable_fp1 = cell(n_c,1);
unstable_fp1 = cell(n_c,1);
stable_fp2 = cell(n_c,1);
unstable_fp2 = cell(n_c,1);

regime = strings(n_c,1);

%% Initial Guesses
range_ticks = [0,1e2,1e4];
range_hosts = [0 40];

L_range = range_ticks;
N_range = range_ticks;
A_range = range_ticks;
LI_range = range_ticks;
NI_range = range_ticks;
AI_range = range_ticks;

IH_range = range_hosts;

%% Start parallel computing
if isempty(gcp('nocreate'))
    parpool;   % start parallel pool only if one doesn't exist
end

unique_solutions1_prev = [];
unique_solutions2_prev = [];

%% Loop over varied param

for ii = 1:n_c
    c = c_values(ii);

    bL = zeros(n_hosts,1);
    bN = zeros(n_hosts,1);
    
    idx = [1:(n_hosts-2) n_hosts];
    
    bL(idx) = b_hat_L(idx)*c;
    bN(idx) = b_hat_N(idx)*c;
    
    bL(n_hosts-1) = b_hat_L(n_hosts-1)*(1-c);
    bN(n_hosts-1) = b_hat_N(n_hosts-1)*(1-c);

    SLsum = 0;
    SNsum = 0;
    
    for i=1:n_hosts-1
        SLsum = SLsum + bL(i)*(H(i)+Lambda(i)/k(i));
        SNsum = SNsum + bN(i)*H(i);
    end
    
    SLsum = SLsum + bL(n_hosts)*(H(n_hosts)*exp(-mu(n_hosts)/2)+Lambda(n_hosts));
    SNsum = SNsum + bN(n_hosts)*H(n_hosts)*exp(-mu(n_hosts)/2);
    
    SAsum = bA*H(n_hosts);

    pL = zeros(n_hosts,1);
    pN = zeros(n_hosts,1);
    
    pL_denom = 0;
    pN_denom = 0;
    
    for i=1:n_hosts-1
        pL_denom = pL_denom + bL(i)*(H(i)+Lambda(i)/k(i))*exp(-mu(i)/2);
        pN_denom = pN_denom + bN(i)*H(i);
    end
    
    pL_denom = pL_denom + bL(n_hosts)*H(n_hosts);
    pN_denom = pN_denom + bN(n_hosts)*H(n_hosts)*exp(-mu(n_hosts)/2);
    
    for i=1:n_hosts-1
        pL(i)=bL(i)*(H(i)+Lambda(i)/k(i))*exp(-mu(i)/2)/pL_denom;
        pN(i)=bN(i)*H(i)/pN_denom;
    end
    
    pL(n_hosts)=bL(n_hosts)*H(n_hosts)/pL_denom;
    pN(n_hosts)=bN(n_hosts)*H(n_hosts)*exp(-mu(n_hosts)/2)/pN_denom;

    params = {r, K, sL, sN, sA, H(1:n_hosts),zL, zN, zA, Lambda, SLsum, SNsum, SAsum, mu, k, betaL, betaN, alpha, pL, pN};

    solutions1 = [];
    solutions2 = [];

    parfor l = 1:length(L_range)
        
        local_sol1 = [];
        local_sol2 = [];

        stopFlag1 = false;
        stopFlag2 = false;

        for n = 1:length(N_range)
        for a = 1:length(A_range)

            L0 = L_range(l);
            N0 = N_range(n);
            A0 = A_range(a);

        for li = 1:length(LI_range)

            if LI_range(li) > L0
                continue
            end

        for ni = 1:length(NI_range)

            if NI_range(ni) > N0
                continue
            end

        for ai = 1:length(AI_range)

            if AI_range(ai) > A0
                continue
            end

        for h1 = 1:length(IH_range)
        for h2 = 1:length(IH_range)
        for h3 = 1:length(IH_range)
        for h4 = 1:length(IH_range)
        for h5 = 1:length(IH_range)

            if IH_range(h1) > H(1), continue, end
            if IH_range(h2) > H(2), continue, end
            if IH_range(h3) > H(3), continue, end
            if IH_range(h4) > H(4), continue, end
            if IH_range(h5) > H(5), continue, end

            fprintf('%d %d %d\n',ii,l,h5)

            x0 = [
                L0
                N0
                A0

                LI_range(li)
                NI_range(ni)
                AI_range(ai)

                IH_range(h1)
                IH_range(h2)
                IH_range(h3)
                IH_range(h4)
                IH_range(h5)
            ];

            [x1,f1,ef1] = fsolve(@(x) mySystemInfection(x,params),x0,opts);
            [x2,f2,ef2] = fsolve(@(x) mySystem2Infection(x,params),x0,opts);
            
            % ----- system 1 -----
            if ~stopFlag1
                if ef1>0 && all(x1>=0) && norm(f1)<1e-6
                    local_sol1 = [local_sol1; x1'];
                end
            end
            
            % ----- system 2 -----
            if ~stopFlag2
                if ef2>0 && all(x2>=0) && norm(f2)<1e-6
                    local_sol2 = [local_sol2; x2'];
                end
            end

        end
        end
        end
        end
        end

        end
        end
        end
        end
        end

        solutions1 = [solutions1; local_sol1];
        solutions2 = [solutions2; local_sol2];

    end

    %% unique solutions
    unique_solutions1 = unique(round(solutions1,3),'rows');
    unique_solutions2 = unique(round(solutions2,3),'rows');

    % If system 1 has solutions, use them as references
    if ~isempty(unique_solutions1)
        ref_fp = unique_solutions1;
    else
        ref_fp = [];                  % no reference solutions
    end
    % Now compare system 2 to system 1
    if ~isempty(unique_solutions2)

        if isempty(ref_fp)
            % Nothing to compare: keep all system 2 solutions
            nonmatching_fp2{ii} = unique_solutions2;
        else
            % Compute differences between each system 2 solution and each system 1 solution
            % Using pdist2 so all pairwise distances are computed
            D = pdist2(unique_solutions2, ref_fp);
            % Keep system 2 points that do NOT match ANY system 1 point
            nonmatch_mask = min(D, [], 2) > 1e2;
            nonmatching_fp2{ii} = unique_solutions2(nonmatch_mask, :);
        end
    else
        nonmatching_fp2{ii} = [];
    end
    %% Stability
    stable1 = [];
    unstable1 = [];
    stable2 = [];
    unstable2 = [];

    % Fixed Point
    for s1 = 1:size(unique_solutions1,1)

        x = unique_solutions1(s1,:)';
        J_fun1 = makeJacobianInfection([params,x(1),x(2),x(3)]);
        J = J_fun1(x);
        eigvals = eig(J);
        if max(abs(eigvals)) < 1
            stable1 = [stable1; x'];
        else
            unstable1 = [unstable1; x'];
        end
    end

    % 2-cycle
    % System 2 stability (ONLY nonmatching points)
    these_fp2 = nonmatching_fp2{ii};       % <- filtered list
    for s = 1:size(these_fp2,1)
        % Use this point as the representative
        x1 = these_fp2(s,:)';
        J_fun1 = makeJacobianInfection([params,x1(1),x1(2),x1(3)]);
        % Generate the orbit in the correct order
        x2 = tickMapInfection(x1,params);
        % Compute Jacobians
        J1 = J_fun1(x1);
        J2 = J_fun1(x2);
        % Jacobian of one trip around the cycle
        M = J2*J1;
        eigvals = eig(M);
        if max(abs(eigvals)) < 1
            stable2 = [stable2; x1'];
        else
            unstable2 = [unstable2; x1'];
        end
    end

    stable_fp1{ii}   = stable1;
    unstable_fp1{ii} = unstable1;
    stable_fp2{ii}   = stable2;
    unstable_fp2{ii} = unstable2;

end

%% Plotting
n = length(c_values);

fp_branch_stable = nan(n,6);
fp_branch_unstable = nan(n,6);

vals_stable = cell(n,1);
vals_unstable = cell(n,1);

lower_branch_unstable_sync = nan(n,6);
upper_branch_unstable_sync = nan(n,6);

for i = 1:n
    
    % ---- Stable FP1 ----
    if ~isempty(stable_fp1{i})
        fp_branch_stable(i,:) = stable_fp1{i}(end,1:6);
    end

    % ---- Unstable FP1 ----
    if ~isempty(unstable_fp1{i})
        for s = 1:size(unstable_fp1{i},1)
            temp_full = unstable_fp1{i}(s,1:6);
            temp = temp_full(1:6);

            if all(temp > 0)
                fp_branch_unstable(i,:) = temp_full;
                break
            end
        end
    end
end
for i = 1:n

    if isempty(unstable_fp2{i})
        continue;
    end
 
    temp_full = unstable_fp2{i}(:,1:6);
    temp_key  = temp_full(:,1:6);

    if isempty(temp_full)
        continue;
    end
    
    [~, ia] = uniquetol(temp_key,1e-2,'ByRows',true);
    vals_unstable{i} = temp_full(ia,:);

    k = size(vals_unstable{i},1);

    if k == 1
        if vals_unstable{i}(1,1) == 0
            lower_branch_unstable_sync(i,:) = vals_unstable{i}(1,:);
        elseif vals_unstable{i}(1,1) > 0
            upper_branch_unstable_sync(i,:) = vals_unstable{i}(1,:);
        end

    elseif k == 2
        vals_unstable{i} = sortrows(vals_unstable{i},1);
        if vals_unstable{i}(2,1) == 0
            lower_branch_unstable_sync(i,:) = vals_unstable{i}(2,:);
            upper_branch_unstable_sync(i,:) = vals_unstable{i}(1,:);
        elseif vals_unstable{i}(1,1) == 0
            lower_branch_unstable_sync(i,:) = vals_unstable{i}(1,:);
            upper_branch_unstable_sync(i,:) = vals_unstable{i}(2,:);
        end
    elseif k == 3
        if vals_unstable{i}(1,1) == 0
            lower_branch_unstable_sync(i,:)  = vals_unstable{i}(2,:);
        end
    elseif k == 4
        lower_branch_unstable_sync(i,:)  = vals_unstable{i}(2,:);
        upper_branch_unstable_sync(i,:)  = vals_unstable{i}(4,:);
    end
end

% Plotting
varNames = {'Larvae (larvae/ha)', 'Nymphs (nymphs/ha)', 'Adults (adults/ha)', 'DIL (larvae/ha)', 'DIN (nymphs/ha)', 'DIA (adults/ha)'};
for col = 1:6        % 1=L, 2=N, 3=A

figure; hold on;
xlabel('$c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel(varNames{col},'FontSize',14);
hold on

% y = lower_branch_unstable_sync(:,col);
% mask = ~isnan(y);
% plot(c_values(mask), y(mask), 'r--','LineWidth',3)
% plot(c_values(mask), y(mask),'w','LineWidth',1)
% 
% y = upper_branch_unstable_sync(:,col);
% mask = ~isnan(y);
% plot(c_values(mask), y(mask),'r--','LineWidth',3)
% plot(c_values(mask), y(mask),'w','LineWidth',1)

plot(c_values, fp_branch_stable(:,col),'b','LineWidth',1.5);
% plot(c_values, fp_branch_unstable(:,col),'r--','LineWidth',1.5);
legend('DIN', 'FontSize', 14, 'Location', 'best');
% legend([hStable, hUnstable], {'Stable', 'Unstable'}, 'FontSize',14,'Location','best');
end

%Plotting NIP

figure; hold on;
xlabel('$c$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('NIP (%)','FontSize',14);
% hUnstable = plot(nan, nan, '--r', 'LineWidth', 1.5);
% hStable = plot(nan, nan, 'b', 'LineWidth', 1.5);
hold on

% y_lower = lower_branch_unstable_sync(:,2);
% y_upper = lower_branch_unstable_sync(:,5);
% mask = ~isnan(y_lower) & ~isnan(y_upper);
% y_avg = y_upper./y_lower;
% plot(c_values(mask), y_avg(mask)*100, 'r--', 'LineWidth', 3)
% plot(c_values(mask), y_avg(mask)*100, 'w', 'LineWidth', 1)

plot(c_values, fp_branch_stable(:,5)./fp_branch_stable(:,2)*100,'b','LineWidth',1.5);
% plot(c_values, fp_branch_unstable(:,col),'r--','LineWidth',1.5);


% legend([hStable, hUnstable], {'Stable', 'Unstable'}, 'FontSize',14,'Location','best');
legend('NIP', 'FontSize',14,'Location', 'best');

delete(gcp('nocreate'));