%% MPCC_path
% solves HJB backwards in time and KFE forward in time
%
% Note:
%   - Relative to the previous version, admits lump-sum transfers

% Initializing Variables
V = V_ss ;  % now call value function V 
mu= mu_ss;
sigma=sigma_ss;

% Emptying Variables
mu_t = zeros(N,T);  sigma_t = zeros(N,T); c_t = zeros(N, T); % just record-keeping   
c    = zeros(N,1);
for tt = T:(-1):1,
    
    % Update current policies
    MPCC_policypath_iv;
    
    % Derivative and Policy Rule
    V1 = (V(2:end) - V(1:end-1))./(s_vec(2:end) - s_vec(1:end-1));
    p=1;
    MPCC_backcons; 
    clear p;
    
    % Dynamics and Payoff Calculations: mu, s, g 
    U =URF(c)                                                         ;
    mu = zeros(N,1)                                                   ;  
    mu(2:end) = (r-pi_t(tt))*s_vec(2:end) + w2 - c(2:end)+tau_rate     ;  
    mu(1) = (r-pi_t(tt))*s_vec(1) + w1 - c(1)+tau_rate                ; 
    sigma = zeros(N,1);  sigma(2:end-1) = s2                          ; 
    
    % record-keeping, could keep track of the value function as well 
    mu_t(:,tt) = mu; sigma_t(:,tt) = sigma; c_t(:,tt) = c;
     
    % Compute Dynamic_HJB
    V=HJB_dynamic_implicit(V,U,s_vec,mu,sigma,rho,dt);
    if plotiter==1
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
    mu = mu_t(:,tt); sigma = sigma_t(:,tt); 
    
    % Solve a forward shift of the KFE 
    [f]=KFE_fs_implicit(f,mu,sigma,s_vec,dt_f,Titer);

    % Record Keeping
    f_t(:,tt) = f;  
    
    if plotiter==1
        figure(pic);
        plot(f); hold on; drawnow;
    end
end
% f_t=[f_0 f_t(:,1:T-1) f_ss];
f_t=[f_0 f_t(:,1:T-1)];
if plotiter==1
    figure(5);
    plot(f_t(:,T-1),'r'); hold on;
    plot(f_ss,'k--','LineWidth',1); hold off; drawnow;
end

%% Check KFE conditions
check=zeros(T,1);
for tt = 1:T-1,
    check(tt)=f_t(:,tt)'*mu_t(:,tt);
end

%% Computing Excess Demand and Supply Functions
% aggregate endowment and consumption 
Y_t = zeros(1,T);  
C_t = zeros(1,T); 
Z_y = zeros(1,T); 
Z_m = zeros(1,T);

% Market clearing along path
for tt = 1:T,
        % Total Endowment people receive, given policy, from distribution
        Y_t(tt) = sum(f_t(:,tt))*w2 + f_t(1,tt)*(w1 - w2);

        % Total Consumption people demand, given policy 
        C_t(tt) = (f_t(:,tt)'*c_t(:,tt));  
        
        % Excess Savings Supply
        Z_y(tt) = (Y_t(tt)-C_t(tt))./Y_t(tt)*100;      
end

% Check Walras law

% Check Steady States
steadycheck=0;
if steadycheck==1
    C_check=(C_t(end)-C_ss)/C_ss;
    M_check=Z_m(end);
end