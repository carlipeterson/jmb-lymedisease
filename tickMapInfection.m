function xnext = tickMapInfection(vars,params)

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

    n_hosts = length(Lambda);

    L = vars(1);
    N = vars(2);
    A = vars(3);
    LI = vars(4);
    NI = vars(5);
    AI = vars(6);
    
    HI = vars(7:6+n_hosts-1);

    tol = 1e-12;

    if A <= tol
        birth_term = 0;
    else
        birth_term = (r*A/2)/((A/2)/H(n_hosts)+K);
    end
    %% -------- SLsum / SNsum / SAsum --------
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

    if L <= tol || SN <= tol
        SA = 1;
    else
        SA = 1-exp(-zA*SAsum/(sN*SN*L));
    end

    % infection dynamics
    infectionLsum = 0; infectionNsum = 0;
    infect_exp = ones(n_hosts-1,1);

    for i = 1:n_hosts-1
        infect_exp(i) = exp(-pN(i)*alpha(i)*sN*SN*LI/H(i));
        infectionLsum = infectionLsum+pL(i)*betaL(i)*(HI(i))/(H(i)+Lambda(i)/k(i));
        infectionNsum = infectionNsum+pN(i)*betaN(i)*HI(i)*exp(-mu(i))/H(i);
    end
    
    xnext = zeros(6,1);

    xnext(1) = birth_term*sL*SL;
    xnext(2) = sN*SN*L;
    xnext(3) = sA*SA*sN*SN*L;
    xnext(4) = birth_term*sL*SL*(infectionLsum);
    xnext(5) = sN*SN*(LI+(L-LI)*(infectionNsum));
    xnext(6) = sA*SA*sN*SN*(LI+(L-LI)*(infectionNsum));

    for i = 1:n_hosts-2
        xnext(i+6) = ((HI(i)+((H(i)-HI(i)+Lambda(i)/2*(1+exp(mu(i)/2)))*(1-infect_exp(i))))*exp(-mu(i)));
    end
    xnext(6+n_hosts-1) = (HI(n_hosts-1)+((H(n_hosts-1)-HI(n_hosts-1))+Lambda(n_hosts-1))*(1-infect_exp(n_hosts-1)))*exp(-mu(n_hosts-1));

end