function [res,aux,D_t]=MPCC_path_solve_FP(rs_t,init_con,steady,paths)
% Declare global variables
MPCC_globals;

% initial condition for distribution 
p_t=rs_t*0+1;
f_0=init_con.f_0;

% Initial Value function
V=steady.V_ss;
f_ss=steady.f_ss; % Used only to draw a picture for the convergence

% Update Borrowing Rate
rb_t=rs_t;
rsp_t = rb_t - rs_t;

RT_t = zeros(1,T)          ;
TE_t = RT_t                ;

% Solve again - with the new price
MPCC_path           ;

% Plot Iteration
Periods = [1:1:T];
if plotit==1
    if rand<0.01
        if rand<0.04
            close all;
        end
        figure(cc);
        subplot(3,1,1);
        plot(Periods, Z_M); 
        title('Excess Money Supply - (Convergence)'); ylabel('Residual over Level'); grid on;
        hold on; drawnow;
        subplot(3,1,2);
        plot(Periods, Z_Y); 
        title('Excess Goods Supply - (Convergence)'); ylabel('Residual over Level'); grid on;
        hold on; drawnow;
        subplot(3,1,3);
        plot(Periods, Z_S); 
        title('Excess Asset Supply - (Convergence)'); ylabel('Residual over Level'); grid on;
        hold on; drawnow;
        figure(cc+1);
        plot(rs_t); drawnow;
    end
end

% Compute Invariant distribution again
if strcmp(clearcond,'Y');
    res=Z_Y;
elseif strcmp(clearcond,'S');
    res=Z_S;
end
aux = Z_S.*D_t;

end
