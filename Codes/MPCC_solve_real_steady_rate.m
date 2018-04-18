% Solving for real rate...

% Solving Steady State Consumption given candidate rate
cond=2*tol;
while cond>tol, % policy iteration with 5 steps 
    % Update Variables to Steady State
    
    % Save this policy
    iter=iter+1;
    c_in=c;
    
    % Update drift and vol
    mu = zeros(N,1)   ;      % Initialize Vector of Drifts
    sigma = zeros(N,1);      % Initial Vector of Volatilities
    
    % Update Drifts and Volatilities
    % Assuming we work with risk-free technology at borrowing constraint | % Adopting Risky Technology Away from borrowing Constraint
    mu(1)     = r*s_vec(1) + (w1 - c(1))+tau          ;
    mu(2:end) = r*s_vec(2:end) + (w2 - c(2:end))+tau  ;
    sigma(2:end-1) = s2                               ;  % volatility is real  
    U =URF(c)                                         ;  % utility
    
    % Solve HJB_implicit
    V=HJB_implicit(U,s_vec,mu,sigma,rho)          ;
    % -> function should be generalized to admit general terminal conditions,
    
    % Updating of Policy Rule at Steady State
    V1 = (V(2:end) - V(1:end-1))./(s_vec(2:end) - s_vec(1:end-1)); % Approximating Value Function
    p=1;
    MPCC_backcons;
    clear p;
    
    % Update condition
    cond=max(c./c_in-1);
end
% Colecting Steady State Objects
c_ss    = c    ;
sigma_ss= sigma;
mu_ss=mu       ;
V_ss=V         ;



% Solving for Steady-State Distributions
[f_ss]=KFE_ss_implicit(mu_ss,sigma_ss,s_vec,dt);
f_ss=f_ss'                  ; 

% Real Aggregate Endowment and Consumption 
Y_ss=sum(f_ss(:)*w2) + f_ss(1)*(w1 - w2);
C_ss=f_ss'*c_ss                         ;
GoodsZ=(Y_ss-C_ss)/Y_ss*100             ;
RS_ss=f_ss'*s_vec                       ; % Real Savings

% Register Initial Conditions
steady.V_ss1     = V_ss       ; 
steady.mu_ss1    = mu_ss      ;
steady.sigma_ss1 = sigma_ss   ;
steady.f_ss1     = f_ss       ;
steady.C_ss1     = C_ss       ;
steady.Y_ss1     = Y_ss       ;
steady.RS_ss1    = RS_ss      ;