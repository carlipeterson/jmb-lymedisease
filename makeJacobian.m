function J_fun = makeJacobian(params)
% makeJacobian  Create a numeric function handle for the Jacobian
%   J_fun = makeJacobian(params)

    % Unpack parameters
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
    n_hosts = length(Lambda);
    n = 3;
    tick_syms = sym('x',[3 1],'real');
    
    
    L = tick_syms(1);
    N = tick_syms(2);
    A = tick_syms(3);
    
    birth_term = ((r*A/2)/((A/2)/H(n_hosts) + K));
    SN = (1-exp(-zN/L*SNsum));
    
    % Define the system
    F1 = birth_term*sL*(1-exp(-zL/birth_term*SLsum));
    F2 = sN*(1-exp(-zN/L*SNsum))*L;
    F3 = sA*(1-exp(-zA/(sN*SN*L)*SAsum))*sN*(1-exp(-zN/L*SNsum))*L;
    
    F = [F1; F2; F3];

    all_vars = [tick_syms];
    
    J = jacobian(F,all_vars);
    
    J_fun_raw = matlabFunction(J,'Vars',{all_vars});
    
    tol = 1e-12;
    
    J_fun = @(x) J_fun_raw(protect_denoms(x,tol));

    function x2 = protect_denoms(x,tol)
    
    x2 = x;
    
    x2(1) = max(x2(1),tol); % L
    x2(2) = max(x2(2),tol); % N
    x2(3) = max(x2(3),tol); % A
    
    end
end
