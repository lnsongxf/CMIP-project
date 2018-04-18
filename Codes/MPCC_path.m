%% MPCC_path
% Update Income path
cond=2*tol;
Y_ss=steady.Y_ss;
y1_out = paths.mu_w_t(1,:)+steady.RT_ss;
y2_out = paths.mu_w_t(2,:)+steady.RT_ss;
y_out=[y1_out; y2_out];

while cond>tol
    % Update Income
    y=y_out;    
    
    % solves HJB backwards in time and KFE forward in time
    ds = s_vec(2)-s_vec(1); 

    % Emptying Variables
    muF_t = zeros(N,T);  muB_t = zeros(N,T); 
    If_t  = zeros(N,T);  Ib_t = zeros(N,T); 
    sigma_t = zeros(N,T); c_t = zeros(N, T); V_t = zeros(N, T); % just record-keeping   
    c     = zeros(N,1);
    dvF   = zeros(N,1); dvB = zeros(N,1);
    muF   = zeros(N,1)   ;      % Initialize Vector of Drifts
    muB   = zeros(N,1)   ;      % Initialize Vector of Drifts
    
    V=steady.V_ss;
    for tt = T:(-1):1
        V_t(:,tt)  = V;

        % Update current price
        p = p_t(tt); 
        y1=y(1,tt);
        y2=y(2,tt); 
        ss = paths.s_bl_index_t(tt); 

        % Update current policies
        MPCC_exopaths;

        % Derivative and Policy Rule
        % Finite difference approximation
        dvF(1:end-1) = (V(2:end) - V(1:end-1))./(s_vec(2:end) - s_vec(1:end-1)); % Forward
        dvF(N)       = (y2 + r_vec(end)*s_vec(end))^(-gamma);
        dvB(2:end)   = (V(2:end) - V(1:end-1))./(s_vec(2:end) - s_vec(1:end-1)); % Forward
        dvB(1)       = (y1 + r_vec(1)*s_vec(1))^(-gamma);

        MPCC_backcons;
        [cF(1:ss)]=min(cF(1:ss),c_bl(1:ss)); % include borrowing constraint limit
        [cB(1:ss)]=min(cB(1:ss),c_bl(1:ss)); % include borrowing constraint limit    
        cF(cF<0)   = 0; cB(cB<0)   = 0;

        % Update Drifts and Volatilities  
        % Assuming we work with risk-free technology at borrowing constraint | % Adopting Risky Technology Away from borrowing Constraint
        muF       = r_vec.*s_vec + (y2 - cF)*p;
        muF(1:ss) = r_vec(1:ss).*s_vec(1:ss) + (y1 - cF(1:ss))*p ;                                                             
        HcF       = URF(cF) + dvF.*muF;      % forward Hamiltonian


        muB       = r_vec.*s_vec + (y2 - cB)*p;
        muB(1:ss) = r_vec(1:ss).*s_vec(1:ss) + (y1 - cB(1:ss))*p ;                                                             
        HcB       = URF(cB) + dvB.*muB;      % backward Hamiltonian

        c0 = r_vec.*s_vec/p + y2; 
        c0(1:ss)      = r_vec(1:ss).*s_vec(1:ss)/p + y1;
    %     [c0(1:ss)]=min(c0(1:ss),c_bl(1:ss)); % include borrowing constraint limit
        dv0 = c0.^(-gamma);    
        H0  = URF(c0);


        % Handling with non-convexities
        Ineither = (1-(muF>0)) .* (1-(muB<0));
        Iunique  = (muB<0).*(1-(muF>0)) + (1-(muB<0)).*(muF>0);
        Iboth    = (muB<=0).*(muF>=0);
        Ib       = Iunique.*(muB<0) + Iboth.*(HcB>HcF);
        If       = Iunique.*(muF>0) + Iboth.*(HcF>=HcB);
        I0       = Ineither;   

        Ib(N) = 1; If(N) = 0; I0(N) = 0;   
        c = cF.*If + cB.*Ib + c0.*I0;

        sigma = ones(N,1)*s1           ; % Initial Vector of Volatilities
        sigma(ss+1:end) = s2*p       ; % volatility is real  

        U =URF(c);

        c_t(:,tt) = c; muF_t(:,tt) = muF; muB_t(:,tt) = muB; sigma_t(:,tt) = sigma;
        If_t(:,tt) = If; Ib_t(:,tt) = Ib; 

        % Compute Dynamic_HJB
        V=HJB_dynamic_implicit(V,U,s_vec,muF,muB,sigma,rho,dt,Ib,If);
        if plotiter==1
            figure(7)
            plot(V); drawnow; hold on;
        end

    end

    %% Solving for the Distribution
    % Initializing Variables
    f=f_0;

    % now solve forward to find the distribution  
    f_t = zeros(N,T);  

    pic=1    ; 
    for tt = 1:T,
        % distribution is computed using this drift and volatility 
        muF = muF_t(:,tt); muB = muB_t(:,tt); sigma = sigma_t(:,tt); 
        If = If_t(:,tt); Ib = Ib_t(:,tt);

        % Solve a forward shift of the KFE 
        [f]=KFE_fs_implicit(f,muF,muB,sigma,s_vec,dt,Ib,If);

        % Record Keeping
        f_t(:,tt) = f;  

        if plotiter==1
            figure(8);
            plot(f); hold on; drawnow;
        end
    end

    % Construct distribution path
    if plotiter==1
        figure(9);
        plot(f_t(:,T-1)); hold on;
        plot(f_ss,'k--','LineWidth',1); hold on; drawnow;
    end

    %% Pre-Shock levels
    f_t = [f_0 f_t(:,1:T-1)];

    % c_t = [steady.c_ss c_t(:,1:T-1)]; 
    % V_t(:,1) = steady.V_ss;
    % muF_t(:,1) = steady.muF_ss;
    % muB_t(:,1) = steady.muB_ss;
    % sigma_t(:,1) = steady.sigma_ss;

    %% Walras's Law Test - Change in aggregate nominal wealth equals change in outside money
    % This condition, and one market clearing implies the other.
    check=zeros(T,1);
    for tt = 1:T,
        check(tt)=f_t(:,tt)'*muF_t(:,tt)*ds;
    end

    % Computing Excess Demand and Supply Functions
    % aggregate endowment and consumption 
    Y_t = zeros(1,T);  
    C_t = zeros(1,T); 
    D_t = zeros(1,T); 
    B_t = zeros(1,T);
    G_t = zeros(1,T); % Vector for government expenditures
    S_t = zeros(1,T); % Vector for FED profits

    % Residuals
    Z_Y = zeros(1,T); 
    Z_S = zeros(1,T);
    Z_M = zeros(1,T);

    % Market clearing along path
    for tt = 1:T,
            w1   = paths.mu_w_t(1,tt); 
            w2   = paths.mu_w_t(2,tt); 
            % Total Endowment people receive, given policy, from distribution
            ss   = paths.s_bl_index_t(tt)                                    ;
            Y_t(tt) = sum(f_t(:,tt)*ds)*w2 + sum(f_t(1:ss,tt)*ds)*(w1 - w2);

            % Total Consumption people demand, given policy 
            C_t(tt) = f_t(:,tt)'*c_t(:,tt)*ds;  

            % Net Asset Position
            S_t(tt) = f_t(:,tt)'*s_vec*ds;

            % Total Savings (Deposits)
            D_t(tt)=f_t(index_o:end,tt)'*s_vec(index_o:end)*ds;

            % Total Reseres --- assuming they are binding
            M0_t(tt)=varrho*D_t(tt);

            if index_o>1
                B_t(tt)=-f_t(1:index_o-1,tt)'*s_vec(1:index_o-1)*ds;
            else
                B_t(tt)=0;
            end

            % Monetary Aggregates
            M1_t(tt) = D_t(tt)              ;
            MM_st(tt) = M1_t(tt)/M0_t(tt)   ;

            % Asset savings condition
            Z_S(tt)=(D_t(tt)+TE_t(tt)-B_t(tt))./D_t(tt);

            % Asset savings condition
            switch lower(mpregime)
                case {'ror'}
                    M0S_t=M0_t;
                    Z_m(tt)=M0_t(tt)-M0S_t(tt);
                case {'ms'}
                    M0S_t=paths.M0_t; % Switch to Monetary Rule
                    Z_m(tt)=M0_t(tt)-M0S_t(tt);
            end

            % FED profits
            G_t(tt)=0; % Government Expenditures

            % Excess Savings Supply
            Z_Y(tt) = (Y_t(tt)-C_t(tt)-G_t(tt))./Y_t(tt);  
    end
    % Update income paths
    y1_out=paths.mu_w_t(1,:)+RT_t(1:T);
    y2_out=paths.mu_w_t(2,:)+RT_t(1:T);
    y_out=[y1_out; y2_out];

    % Find the optimality condition
    cond=max(max(abs(y-y_out)./y));
end
