%% MPCC_path_solve
% loop over prices
function [res_t]=MPCC_path_solve_iv(pi_t,init_con,steady,paths);
global rho gamma c_bl c_dl;
% global rho gamma iter r c_bl;
global URF Upinv T iter Titer dt dt_f s_vec N pic plotit plotiter;

% Initial Conditions
f_0     = init_con.f_0             ;  % initial conditions

% Variables from earlier computation
V_ss    = steady.V_ss        ; 
mu_ss   = steady.mu_ss       ;
sigma_ss= steady.sigma_ss    ;
s       = s_vec              ;
C_ss     = steady.C_ss       ;
Y_ss     = steady.Y_ss       ;
f_ss     = steady.f_ss       ;
pi_ss    = steady.pi_ss      ;

%% Solution along a transition given guess for prices
mod=pi_t;
arp=1   ;
if arp==1
    
end
pi_t=pi_t;
MPCC_path_rb;

%% Plots
Periods = [1:1:T];
if plotiter==1
    figure(4);
    plot(Periods, Y_t,'b');

    hold on
    plot(Periods, C_t,'r');
end

if plotit==1&&rand<0.01
    figure(10);
    subplot(2,1,1);
    plot(Periods, Z_m,'b'); 
    title('Excess Money Supply - (Convergence)'); ylabel('Residual over Level'); grid on;
    hold on; drawnow;
    subplot(2,1,2);
    plot(Periods, Z_y); 
    title('Excess Goods Supply - (Convergence)'); ylabel('Residual over Level'); grid on;
    hold on; drawnow;
    
    figure(8);
    plot(Periods, pi_t); title('Inflation Rate'); 
end

% Register Residuals
res_t=[Z_y(1:end)];

% res_t=[Z_y(1) Z_m(2:end)];
% res_t=[Z_m(1:2:end-1)+Z_y(2:2:end)];
% res_t=[Z_m(1:end).^2+Z_y(1:end).^2];