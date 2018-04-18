cond=2*tol;
iter = 0;
dvF = zeros(N,1); dvB = zeros(N,1);

while cond>tol % value function iteration
    V_in =V;
    % Finite difference approximation
    dvF(1:end-1) = (V(2:end) - V(1:end-1))./(s_vec(2:end) - s_vec(1:end-1)); % Forward
    dvF(N)       = (y2 + r_vec(end)*s_vec(end))^(-gamma);
    dvB(2:end)   = (V(2:end) - V(1:end-1))./(s_vec(2:end) - s_vec(1:end-1)); % Backward
    dvB(1)       = (y1 + r_vec(1)*s_vec(1))^(-gamma);
    
    MPCC_backcons;
    [cF(s_bl_index)]=min(cF(s_bl_index),c_bl(s_bl_index)); % include borrowing constraint limit
    [cB(s_bl_index)]=min(cB(s_bl_index),c_bl(s_bl_index)); % include borrowing constraint limit    
    cF(cF<0)   = 0; cB(cB<0)   = 0;
    
    % Update drift and vol
    muF = zeros(N,1)   ;      % Initialize Vector of Drifts
    muB = zeros(N,1)   ;      % Initialize Vector of Drifts
    sigma = zeros(N,1);      % Initial Vector of Volatilities
    index=logical(1-s_bl_index);
    
    % Update Drifts and Volatilities  
    % Assuming we work with risk-free technology at borrowing constraint | % Adopting Risky Technology Away from borrowing Constraint
    muF(s_bl_index)     = r_vec(s_bl_index).*s_vec(s_bl_index) + (y1 - cF(s_bl_index))*p ;                                                             
    muF(index) = r_vec(index).*s_vec(index) + (y2 - cF(index))*p  ;
    HcF       = URF(cF) + dvF.*muF;      % forward Hamiltonian
    
    muB(s_bl_index)     = r_vec(s_bl_index).*s_vec(s_bl_index) + (y1 - cB(s_bl_index))*p ;                                                             
    muB(index) = r_vec(index).*s_vec(index) + (y2 - cB(index))*p  ;
    HcB       = URF(cB) + dvB.*muB;      % backward Hamiltonian
   
    c0 = r_vec.*s_vec/p + y2; 
    c0(s_bl_index)      = r_vec(s_bl_index).*s_vec(s_bl_index)/p + y1;
    [c0(s_bl_index)]=min(c0(s_bl_index),c_bl(s_bl_index)); % include borrowing constraint limit
    dv0 = c0.^(-gamma);

    H0  = URF(c0);

    % Handling with non-convexities
    Ineither = (1-(muF>0)) .* (1-(muB<0));
    Iunique  = (muB<0).*(1-(muF>0)) + (1-(muB<0)).*(muF>0);
    Iboth    = (muB<0).*(muF>0);
    Ib       = Iunique.*(muB<0) + Iboth.*(HcB>=HcF);
    If       = Iunique.*(muF>0) + Iboth.*(HcF>=HcB);
    I0       = Ineither;   
    Ib(N) = 1; If(N) = 0; I0(N) = 0;        
    
    c = cF.*If + cB.*Ib + c0.*I0;
    
    sigma(2:end) = s2*p                           ;  % volatility is real      
    U =URF(c);
    
    % Solve HJB_implicit
    V=HJB_implicit(U,s_vec,muF,muB,sigma,rho,Delta,V,Ib,If);
    
    
    % Update condition
    cond=max(abs(V./V_in-1));
    
    % Plot if needed
    if plotiter==1
        plot(V); drawnow; hold on;
    end
    iter=iter+1;
end
% Colecting Steady State Objects
c_ss     = c    ;
sigma_ss = sigma;
muF_ss   = muF   ;
muB_ss   = muB   ;
V_ss     = V    ;

% Reporting Steady State Solutions
% disp('**** Steady State Solution ***');
% disp('Computation time:');
% toc
% disp(['Iterations:   ' num2str(iter)]);
% toc
