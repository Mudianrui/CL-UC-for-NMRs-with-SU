function [sys,x0,str,ts] = xyth_Scontroller(t,x,u,flag)
%% Scontroller
    switch flag
    case 0
        [sys,x0,str,ts]=mdlInitializeSizes;
    case 1
        sys=mdlDerivatives(t,x,u);
    case 3
        sys=mdlOutputs(t,x,u);
    case {2,4,9}
        sys=[];
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
    end


function [sys,x0,str,ts]=mdlInitializeSizes
%% mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 6;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 16;
    sizes.NumInputs      = 25;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0  = [0,0,0,0,0,0];
    str = [];
    ts  = [];
    
    global T
    T = 3;
    
    global alpha td Te cw
    alpha = 10;
    td = 0.5;
    Te = 1.5;
    cw = 1;
    
    blocks = find_system(bdroot,'BlockType','TransportDelay');
    set_param(blocks{1},'DelayTime', num2str(td));
    
    global isReachTe1 isReachTe2 DeltasHat
    isReachTe1 = 0;
    isReachTe2 = 0;
    DeltasHat = [0;0];


function sys=mdlDerivatives(~,sx,u)
%% mdlDerivatives
    t = u(1);
    xs = u(7);
    ys = u(8);
    th = u(9);
    v = u(10);

    global alpha
    E=exp(alpha*t);
    sys(1) = v*cos(th);
    sys(2) = v*sin(th);
    sys(3) = sin(th)*E;
    sys(4) = cos(th)*E;
    sys(5) = (xs-sx(1))*E;
    sys(6) = (ys-sx(2))*E;


function sys=mdlOutputs(~,sx,u)
    t = u(1);
    global controllerTypeK modelTypeAll
    global lc T
    global alpha td Te cw
    global isReachTe1 isReachTe2 ThetaDeltasTe ThetaTe
    global DeltasHat
    global dasz2
    % input
    xr = u(2);
    yr = u(3);
    pr = [xr;yr];
    thr = u(4);
    dxr = u(5);
    dyr = u(6);
    dpr = [dxr;dyr];
    % measurement value
    xs = u(7);
    ys = u(8);
    th = u(9);
    v = u(10);
    w = u(11);
    s = [v;w];
    % intermediate variable
    Delta_s_hatx = u(12);
    Delta_s_haty = u(13);
    Delta_s_hat = [Delta_s_hatx;Delta_s_haty];

    %% Set expected values
    Psi = [cos(th),-sin(th);sin(th),cos(th)];
    Psic = Psi*[1,0;0,lc];
    Psis = Psi;
    dPsis = w*[-sin(th),-cos(th);cos(th),-sin(th)];
    ps = [xs;ys];
    Delta_c = [lc*cos(th);lc*sin(th)];
    
    %% Composite learning estimates the true value of sampling error: calculating prediction error
    switch controllerTypeK
        case {54,56}
            En=exp(-alpha*t);
            Idpo = [sx(1);sx(2)];
            IPsisE = [sx(4),-sx(3);sx(3),sx(4)];
            Ips_poE = [sx(5);sx(6)];
            Psisf = alpha*Psis-alpha^2*En*IPsisE;
            pf = alpha*(ps-Idpo)-alpha^2*En*Ips_poE;
            if t<Te
                dThetaDeltas = Psisf'*pf;
            else
                dThetaDeltas = [0;0];
            end
            ThetaDeltas_te = [u(14);u(15)];
            ThetaDeltas_t = [u(20);u(21)];
            if t<td
                ThetaDeltas = ThetaDeltas_t;
            elseif t<Te
                ThetaDeltas = ThetaDeltas_t-ThetaDeltas_te;
            else
                if isReachTe1==1
                    ThetaDeltas = ThetaDeltasTe;
                else
                    isReachTe1 = 1;
                    ThetaDeltas = ThetaDeltas_t-ThetaDeltas_te;
                    ThetaDeltasTe = ThetaDeltas;
                end
            end
            if t<Te
                dTheta = Psisf'*Psisf;
            else
                dTheta = [0,0;0,0];
            end
            Theta_te = [u(16),u(17);u(18),u(19)];
            Theta_t = [u(22),u(23);u(24),u(25)];
            if t<td
                Theta = Theta_t;
            elseif t<Te
                Theta = Theta_t-Theta_te;
                DeltasHat = pinv(Psisf)*pf;
            else
                if isReachTe2==1
                    Theta = ThetaTe;
                else
                    isReachTe2 = 1;
                    Theta = Theta_t-Theta_te;
                    ThetaTe = Theta;
                end
            end
            epslion = ThetaDeltas-Theta*Delta_s_hat;
            trueDelta = Delta_s_hat+pinv(Theta)*epslion;
    end

    %% Calculate/estimate the coordinates and errors of control points
    switch controllerTypeK
        case {56}
            trnsType = 0;
            trnst = trns(t,Te,T,trnsType);
            dtrnst = trns_derivative(t,Te,T,trnsType);
            pe1 = ps+Delta_c-pr;
            pe2 = ps-Psis*trueDelta+Delta_c-pr;
            pes = (1-trnst)*pe1+trnst*pe2;
            pehat = pes;
        otherwise
            pehat = ps-Psis*Delta_s_hat+Delta_c-pr;
    end
    ex = pehat(1);
    ey = pehat(2);

    %% Error transformation
    switch controllerTypeK
        case {52,54}
            epsilon_ex = 0.1;
            epsilon_ey = 0.1;
            xi_ex_p = 1;
            rou_ex = PPC_rho_inf(t,T,epsilon_ex,xi_ex_p,0);
            xi_ey_p = 1;
            rou_ey = PPC_rho_inf(t,T,epsilon_ey,xi_ey_p,0);
        case {56}
            epsilon_ex = 0.2;
            epsilon_ey = 0.2;
            alphaR = 1;
            amplifierType = 4;
            xi_ex_p = 1;
            delta_ex = 0.1;
            rou_ex = PPC_rho_inf(t,T,epsilon_ex,xi_ex_p,0);
            drou_ex = PPC_rho_inf_derivative(t,T,epsilon_ex,xi_ex_p,0);
            alpha_ex = DirectPPC_amplifier(ex,rou_ex,delta_ex,alphaR,xi_ex_p,amplifierType);
            xi_ex = alpha_ex*ex;
            [mu_ex,upsilon_ex] = DirectPPC_amplifier_derivative(ex,rou_ex,drou_ex,delta_ex,alphaR,xi_ex_p,amplifierType);
            xi_ey_p = 1;
            delta_ey = 0.1;
            rou_ey = PPC_rho_inf(t,T,epsilon_ey,xi_ey_p,0);
            drou_ey = PPC_rho_inf_derivative(t,T,epsilon_ey,xi_ey_p,0);
            alpha_ey = DirectPPC_amplifier(ey,rou_ey,delta_ey,alphaR,xi_ey_p,amplifierType);
            xi_ey = alpha_ey*ey;
            [mu_ey,upsilon_ey] = DirectPPC_amplifier_derivative(ey,rou_ey,drou_ey,delta_ey,alphaR,xi_ey_p,amplifierType);
    end

    %% controller
    switch floor(controllerTypeK/10)
        case 5
            switch mod(controllerTypeK, 10)
                case 2
                    kp = 10;
                    uvw = Psic\(dpr-kp*(ps-Psis*Delta_s_hat+Delta_c-pr));
                    uv = uvw(1);
                    uw = uvw(2);
                    kd = 100;
                    if modelTypeAll==1
                        dDelta_s_hat = -Psis'\(Psic*s-dpr)-kd*Delta_s_hat;
                    else
                        dDelta_s_hat = -Psis'\(Psic*uvw-dpr)-kd*Delta_s_hat;
                    end
                    dDelta_s_hatx = dDelta_s_hat(1);
                    dDelta_s_haty = dDelta_s_hat(2);
                case 4
                    kp = 1;
                    uvw = Psic\(dpr-kp*(ps-Psis*Delta_s_hat+Delta_c-pr));
                    uv = uvw(1);
                    uw = uvw(2);
                    kD = 10;
                    if modelTypeAll==1
                        dot = -Psis'\(Psic*s-dpr)+kD*pinv(Theta)*epslion;
                    else
                        dot = -Psis'\(Psic*uvw-dpr)+kD*pinv(Theta)*epslion;
                    end
                    if norm(Delta_s_hat)<cw
                        Projection = dot;
                    else
                        Projection = dot-Delta_s_hat*Delta_s_hat'/(Delta_s_hat'*Delta_s_hat)*dot;
                    end
                    Gamma = 1*eye(2,2);
                    dDelta_s_hat = Gamma*Projection;
                    dDelta_s_hatx = dDelta_s_hat(1);
                    dDelta_s_haty = dDelta_s_hat(2);
                case 6
                    kp = 5;
                    xi_p = [xi_ex;xi_ey];
                    upsilon_p = [upsilon_ex;upsilon_ey];
                    invmu_p = [1/mu_ex,0;0,1/mu_ey];
                    mu_p = [mu_ex,0;0,mu_ey];
                    Psics_hat = [cos(th),-lc*sin(th)+(1-trnst)*(-sin(th)*Delta_s_hatx-cos(th)*Delta_s_haty);
                                sin(th),lc*cos(th)+(1-trnst)*(cos(th)*Delta_s_hatx-sin(th)*Delta_s_haty)];
                    uvw = Psics_hat\(dtrnst*Psis*Delta_s_hat+dpr-invmu_p*upsilon_p-kp*invmu_p*xi_p);
                    uv = uvw(1);
                    uw = uvw(2);
                    dasz2 = Psic'*mu_p'*xi_p;
                    kD = 1;
                    star = ((1-trnst)*dPsis'-dtrnst*Psis')*mu_p'*xi_p;
                    dot = star+kD*pinv(Theta)*epslion;
                    if norm(Delta_s_hat)<cw
                        Projection = dot;
                    else
                        Projection = dot-Delta_s_hat*Delta_s_hat'/(Delta_s_hat'*Delta_s_hat)*dot;
                    end
                    Gamma = 1*eye(2,2);
                    dDelta_s_hat = Gamma*Projection;
                    dDelta_s_hatx = dDelta_s_hat(1);
                    dDelta_s_haty = dDelta_s_hat(2);
            end
    end

    sys(1) = uv;
    sys(2) = uw;
    sys(3) = thr;
    sys(4) = rou_ex;
    sys(5) = rou_ey;
    sys(6) = 0;
    sys(7) = ex;
    sys(8) = ey;

    switch controllerTypeK
        case {52}
            sys(9) = dDelta_s_hatx;
            sys(10) = dDelta_s_haty;
            sys(11) = 0;
            sys(12) = 0;
            sys(13) = 0;
            sys(14) = 0;
            sys(15) = 0;
            sys(16) = 0;
        case {54,56}
            sys(9) = dDelta_s_hatx;
            sys(10) = dDelta_s_haty;
            sys(11) = dThetaDeltas(1);
            sys(12) = dThetaDeltas(2);
            sys(13) = dTheta(1,1);
            sys(14) = dTheta(1,2);
            sys(15) = dTheta(2,1);
            sys(16) = dTheta(2,2);
    end
