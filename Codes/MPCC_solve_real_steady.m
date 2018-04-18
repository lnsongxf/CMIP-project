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

% Reporting Steady State Solutions
disp('**** Steady State Solution ***');
disp('Computation time:');
toc
disp(['Iterations:   ' num2str(iter)]);
toc
