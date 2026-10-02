function J_fun = makeJacobianInfection(params)

r = params{1};
K = params{2};
sL = params{3};
sN = params{4};
sA = params{5};
H = params{6};
zL = params{7};
zN = params{8};
zA = params{9};
Lambda = params{10};
SLsum = params{11};
SNsum = params{12};
SAsum = params{13};
mu = params{14};
k = params{15};
betaL = params{16};
betaN = params{17};
alpha = params{18};
pL = params{19};
pN = params{20};
L = params{21};
N = params{22};
A = params{23};

n_hosts = length(Lambda);

tick_syms = sym('x',[1 1],'real');
HI = sym('HI',[n_hosts-1 1],'real');

LI = tick_syms(1);


%% -------- birth_term --------

birth_term = ((r*A/2)/((A/2)/H(n_hosts)+K));

%% -------- SLsum / SNsum / SAsum --------
tol = 1e-12;

    if A <= tol
        birth_term = 0;
    else
        birth_term = (r*A/2)/((A/2)/H(n_hosts)+K);
    end
    if birth_term <= tol
        SL = 1;
    else
        SL = 1-exp(-zL*SLsum/birth_term);
    end

    if L <= tol
        SN = 1;
    else
        SN = 1-exp(-zN*SNsum/L);
    end

    if L<=tol || SN<=tol
        SA = 1;
    else
        SA = 1-exp(-zA*SAsum/(sN*SN*L));
    end

%% infection

infectionLsum = 0;
infectionNsum = 0;
infect_exp = sym(ones(n_hosts-1,1));

for i=1:n_hosts-1
    infect_exp(i)=exp(-pN(i)*alpha(i)*sN*SN*LI/H(i));
    infectionLsum = infectionLsum+pL(i)*betaL(i)*(HI(i))/(H(i)+Lambda(i)/k(i));
    infectionNsum = infectionNsum+pN(i)*betaN(i)*HI(i)*exp(-mu(i))/H(i);
end

%% ticks

F1 = ((r*A/2)/((A/2)/H(n_hosts)+K))*sL*SL*infectionLsum;

F_ticks = [F1];

%% hosts

F_hosts = sym(zeros(n_hosts-1,1));

for i=1:n_hosts-2
    F_hosts(i) = ((HI(i)+((H(i)-HI(i)+Lambda(i)/2*(1+exp(mu(i)/2)))*(1-infect_exp(i))))*exp(-mu(i)));
end
F_hosts(n_hosts-1) = (HI(n_hosts-1)+((H(n_hosts-1)-HI(n_hosts-1))+Lambda(n_hosts-1))*(1-infect_exp(n_hosts-1)))*exp(-mu(n_hosts-1));

F = [F_ticks;F_hosts];
all_vars = [tick_syms;HI];

infected_vars = [LI;HI];
F_infected = [F1;F_hosts];

J = jacobian(F_infected,infected_vars);

J_fun = matlabFunction(J,'Vars',{infected_vars});

end
