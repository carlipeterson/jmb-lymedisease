function xnext = tickMap(vars,params)

    L = vars(1);
    N = vars(2);
    A = vars(3);

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

    birth_term = (r*A/2)/((A/2)/H(n_hosts)+K);
    
    xnext = zeros(3,1);

    xnext(1) = birth_term*sL*SL;
    xnext(2) = sN*SN*L;
    xnext(3) = sA*SA*sN*SN*L;

end