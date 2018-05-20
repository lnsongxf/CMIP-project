%% Moll's algorithm
% 
cond =2*tol_dyn;
r_t = rs_t_o;
iter=0;
dZ_S = nan(T,1);
% xi   = 1e-5;
D = [];
r = [];
Z = [];
% flagx10 =0;

while cond>tol_dyn
    [~,Z_S,D_t] = eqpath(r_t);

    r_t        = r_t - xi.*Z_S;  

    
    if iter<200
%         r_t(T-10:T)=rs;                 % We know that r will approach the stationary value, so we will force smoothing around the stationary value.
        r_t(10:T)=smooth(r_t(10:T),20); % Smooth updated r
    elseif iter<300
%         r_t(T-5:T)=rs;
        r_t(5:T)=smooth(r_t(5:T),10);
    elseif iter<400
%         r_t(T-5:T)=rs;
        r_t(5:T)=smooth(r_t(5:T),5);
%     elseif max(abs(Z_S))<1/2 & xi==1e-5
%         xi=xi*10;        
    end   
%     r_t(T-5:T)=rs;
    
%     if max(abs(Z_S))<1 && flagx10 == 0;
%         xi = xi*10;
%         flagx10=1;
%     end
    
    if mod(iter,10)==0
        figure(cc);        
        subplot(3,1,1);
        plot(r_t);   title('r_t');   
        subplot(3,1,2);
        plot(Z_S);  title('D_t-B_t');
        subplot(3,1,3);
        plot(Z_S./D_t);  title('res_t');
        drawnow;
        if mod(iter,50), clf; end
    end
    
    cond = max(abs(Z_S./D_t));
    iter=iter+1;
    D(:,iter) = D_t;
    r(:,iter) = r_t;
    Z(:,iter) = Z_S;
end

