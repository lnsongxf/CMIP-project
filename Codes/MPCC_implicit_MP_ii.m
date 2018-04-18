%% MPCC_implicit

% Main Code to Solve and Simulate the MPPC model
% (c) Bigio and Sannikov
clear; 
close all;

%% Notes
% - Relative to previous version (version iv):
%     - Includes Monetary Policy with FED profits thrown to the sea
%     (instead of lump-sum)
%     - Monetary regime can either fix rate, or fix money supply
%     - Rate on reserves pins down borrowing-lending spread
% 
% - Missing things we spoke about:
%     - from this code...
%     - Includes technology choice
%     - interest rate - add arbitrage to market clearing check
%     - bliss point in consumption - make it work for b
%     - Generalization of aggregate shocks - event trees
%     - Lump-sum monetary transfers -> add condition
%     - Default
%     - Output feedback loop?

%% Globals
MPCC_globals_MP; % Define list of globals

%% Code Parameters
% [-] Move solution Method choices here
% Options for Solution of ODE ---> Used for HJB
plotiter=0; % set to 1, if want to see algorithm converge
plotit  =0; % 
breakit =0; % 
printit =0; % 

% Plot Preferences
ser1Width=1;

%% Code Parameters
% Approximate Amount of grid points
N = 150        ; % Number of Gridpoints in Space

% Time Parameters
T     = 60    ; % Vector of Intervals
Titer = 1     ; % Iterations

% From older Code
T_pre = 20       ; % Time Lapse prior to transition
T_post= 40       ; % Time Lapse after transition
Tspan = [1/12 0] ; % Vector Time, within interval

% State-Space Grid
dt = 1    ; % Time Dimention - Monthly Model (12)
dt_f = dt ; % A step size for the KFE

% Optimization Parameters 
tol=10e-6;
cond=2*tol;
options = optimset('TolFun',1e-4);

% Plot Properties
s_plotmax=N/2;
t_plotmax=60 ;

color1_ss   =[0.8 0.8 0.8];
color1_shock=[0.8 0.2 0.2];
color2_shock=[0.6 0.2 0.2];

%% Model Parameter Definitions
% [-] Pass Labels to Functions

% Auxiliaries
iter=0;

% Model Parameters - At Steady State
rho = 0.1   ; % agent's individual discount rates -> presumably need rho >= rs to make sure savings don't blow up
gamma = 0.4 ; % agent's risk aversion
w1 = 2      ; % mean low process
w2 = 10     ; % mean high process 
s1 = 0      ; % volatility low process
s2 = 5      ; % volatility 

% Borrowing Limit - nominal
s_bar  = -10*w1  ;  % debt limit
s_bl   = s_bar   ;  % borrowing limit - below this point, you can accumulate interest but not borrow

% Vectorizations
mu_w=[w1 w2]'    ;
sigma_w=[s1 s2]' ;

% Bliss point
c_bliss=-w1*0; % Potentiall bliss point to break homogeneity assumption - useful to undo risk-aversion effect

%% Functional Forms 
if gamma~=1
    URF=@(c) ((c-c_bliss).^(1-gamma)-1)/(1-gamma); % Utility Return Function
    Upinv=@(c) (c-c_bliss).^(-1/gamma); % Utility Return Function
else
    URF=@(c)   log(c-c_bliss); % Utility Return Function
    Upinv=@(c) exp(c-c_bliss); % Utility Return Function
end

%% Policy and Assumed Prices
mpregime='ror' ; % target overnight rate. set liq to 0, to be in regime with rate on excess reserves only
% set to:
% -------------------------------------------------------------------------
%
% ror: If based on payments on interest rates on reserves, fed accomodates
% to policy
%
% ms: money supply follows an exogenous path holding corridor rates fixed exgogenous
%
% -------------------------------------------------------------------------
resreq=0.1       ; % Reserve Requirements
liq   =resreq*2  ; % Payments Need - Average discount window loans over reserves
rr    = 0.00     ; % Rate on Excess Reserves
rsp   = 0.05     ; % Width of Corridor Window
p    = 1         ; % Initial Steady-State Price Level 

% lending rate from zero-profit condition - define inside whenever called
r_b_zp=@(rl,rr,rsp) (rl-resreq*rr+liq*(rsp+rr))/(1-resreq);

%% Deduced Parameters  
s_max = 1/2*w2/rho                      ;  % consider savings up to 4 times PV of expected risky income 

%% Discretization of State Space
% Savings space and consumption space Gus
s_vec  = linspace(s_bar, s_max,N)'           ; % Space of (nominal) savings 

% Initialized Variables
mu_ss  = zeros(1,N); % Drift at Steady State
sigma_ss = zeros(1,N); % Volatility at Steady State

% Borrowing Limit - Static ---- NEED TO ADAPT HERE, to admit credit crunch
% shock
c_dl = rho*s_bar + w1; % Stays Constant - 

% Borrowing Limit - Dynamic
c_bl(2:N,1)= w2  ;

%% Experiments - Determinisitic paths
% Only allowing Technology / policies to change 
% [1.a ] *Means and Volatilities along Transition
shock_mu   =1.0                             ;
mu_w_t     =mu_w*ones(1,T)                  ;
mu_w_t(2,T_pre+1:T_post)= mu_w(2)*shock_mu  ;

% [1.b ]Means and *Volatilities along Transition
sigma_w_t=sigma_w*ones(1,T)  ;

% path
shock_sigma=1.0;
sigma_w_t(2,T_pre+1:T_post)= s2*shock_sigma;

% [2.a] Vector of Policy Rates - interest on reserves
rr_shock                     = rr    ;
rr_t(1:T)                    = rr   ; % Interest on reserves
rr_t(T_pre+1:T_post)=rr_shock       ;

% [2.a] Vector of Policy Rates - interest on reserves
rsp_shock                     = rsp*2  ;
rsp_t(1:T)                    = rsp    ; % Interest on reserves
rsp_t(T_pre+1:T_post)=rsp_shock        ;

% [3] Borrowing Limits 
shock_s_bl=1.0  ;
s_bl_t(1:T)=s_bl;
s_bl_t(T_pre+1:T_post)=s_bl*shock_s_bl;

% [4] Transfer Shocks - Expected
shock_T=0.00                        ; % Rate of money growth do to transfers 
T_t(1:T)=1                          ;
T_t(T_pre+1:T_post)=(1+shock_T*dt_f);

% [5] Helicopter Drop
Delta_ingrid=1                      ; % Unexpected shock

% write down paths
paths.rr_t=rr_t          ;
paths.rsp_t=rsp_t         ;
paths.mu_w_t=mu_w_t      ;
paths.sigma_w_t=sigma_w_t;
paths.s_bl_t=s_bl_t      ;
paths.T_t=T_t            ; 

% Paths for exogenous 
paths_list={'rr_t','rsp_t','mu_w_t','sigma_w_t','s_bl_t','T_t'};
paths_list_tit={'r^o^r_t','r^d^w_t','\mu_2','\sigma_2','s^b^l_t','T_t'};

%% Add Separate code to describe shocks
MPPC_report_shocks

%% [1] Solving Initial Steady State
% Set timers
tic
iter=0;

% Initial price is p=1
p=1                     ;
p_ss=p                  ;
rs_o=0.04               ;

% Solve using unit price
inputs.p=1; % Normalization of price level

% Initial conjecture for consumption
c_vec = c_bliss+0.8*linspace(w1, w2, N)'+rho*s_vec; % Initial guess - permanent income hypothesis
c_vec = max(c_vec,0.1)                           ; % 

%Solution of HJB Equation using Policy Iteration Method
c=c_vec;
inputs.c=c;

% res = MPCC_solve_steady_MP(rs_o,targets,paths);
plot_zA=0;
if plot_zA==1
    fplot(@(x) MPCC_solve_steady_MP(x,inputs,paths),[-0.01 0.05]); grid on;
end

% Define and solve steady state rates
res=@(rs_o) MPCC_solve_steady_MP(rs_o,inputs,paths);
rs_ss=fsolve(res,rs_o,options);

% Update Rate
% rs=exp(rs)                ; % allows only positive rates     
rs=rs_ss                    ; % impose borrowing limit
rb=r_b_zp(rs,rr,rsp)        ; % lending rate
rb_ss=rb                    ; % lending rate - steady state
rr_ss=rr                    ; % reserve rate - steady state

%-----------------------------------
% Recomputing solution at steady state
%-----------------------------------
% Construct Interest Rate Vector
r_vec=s_vec       ; % Initialize it
r_vec(s_vec>0)=rs ;
r_vec(s_vec<=0)=rb;

% Call Solver
MPCC_solve_steady;

% Reporting Steady State Solutions
disp('**** Steady State Solution ***');
disp('Computation time:');
toc
disp(['Iterations:   ' num2str(iter)]);
toc

% [2] Call Solve Values targets
% Solving for Steady-State Distributions
[f_ss]=KFE_ss_implicit(mu_ss,sigma_ss,s_vec,dt) ;
f_ss=f_ss'                                      ; 

% Aggregate endowment and consumption 
Y_ss=sum(f_ss(:)*w2) + f_ss(1)*(w1 - w2)        ;
C_ss=f_ss'*c_ss                                 ;

% Savings+Borrowings=0
S_ss=f_ss'*s_vec                                ; % Nominal Money Stock

% Compute Deposits
index_o=find(s_vec>=0,1,'first')                ;

% Computing Savings and Borrowings
D_ss=f_ss(index_o:end)'*s_vec(index_o:end)      ;
M0_ss=resreq*D_ss                               ;
if index_o>1
    B_ss=-f_ss(1:index_o-1)'*s_vec(1:index_o-1) ;
else
    B_ss=0                                      ;
end

% Monetary Aggregates
M1_ss = D_ss        ;
MM_ss = M1_ss/M0_ss ;

if plotiter==1
    plot(s_vec,f_ss); xlim([s_vec(1) s_vec(end)]); ylim([0 max(f_ss)]); drawnow;
end
AssetsZ=D_ss-B_ss;

% FED profits
FED_profs=(rb-rr)*M0_ss+(rsp+rr)*liq*D_ss;

% Including FED profits
% Goods Market Clearing - Doesn't hold because FED makes profits GoodsZ=(Y_ss-C_ss)/Y_ss            ;
GoodsZ=(Y_ss-C_ss-FED_profs)/Y_ss ;

%% Saving steady state variables and initial conditions
% Recording Terminal Steady State
steady.V_ss      = V_ss      ; 
steady.mu_ss     = mu_ss     ;
steady.sigma_ss  = sigma_ss  ;
steady.p_ss      = p_ss      ;
steady.MStock_ss = M1_ss     ;
steady.MStock_ss = M0_ss     ;
steady.MStock_ss = MM_ss     ;
steady.f_ss      = f_ss      ;
steady.C_ss      = C_ss      ;
steady.Y_ss      = Y_ss      ;

% Initial conditions - distribution
f_0=f_ss;
init_con.f_0=f_0;

% If policy regime is money rule - Adapt here if money path changes
switch lower(mpregime)
    case {'ror'}
        
    case {'ms'}
        paths.M0_t=M0_ss*ones(1,T);
end

%% Solving Transition 
% Guess for prices
p_t_o  = ones(1, T)    ;
rs_t_o = rs_ss*ones(1, T)   ;
q_vec_o= [p_t_o; rs_t_o]; % Guess of p_t and r_d in vectorial form

% Solving for the price path
res_sup=@(q_vec) MPCC_path_solve_MP_ii(q_vec,init_con,steady,paths);

% Solving the main function
tic
q_t_sol=fsolve(res_sup,q_vec_o,options);
disp('Time to Find Solution');
toc

% Calling function to compute objects
p_t  = q_t_sol(1,:); % Update final price solution
rs_t = q_t_sol(2,:); % Update final rate solution

% Save variables
rb_t =r_b_zp(rs_t,rr_t,rsp_t); % use formula

% Initial Value function
V=V_ss;

% Solve again - with the new price
MPCC_path_MP           ;
return
%% Transition Plots
% s_plotmax - max s for plots
% t_plotmax - max s for plots
% color1_ss   =[0.8 0.8 0.8];
% color1_shock=[0.8 0.2 0.2];
% color2_shock=[0.6 0.2 0.2];


% Compute Solution
Periods = (1:1:T);

% Plotting Residual Functions
figure;
plot(Periods, Z_m,'b'); title('Excess Money Supply'); drawnow;

figure;
plot(Periods, Z_y,'b'); title('Excess Supply of Goods'); drawnow;

figure;
plot(Periods, Z_a,'b'); title('Excess Supply of Assets'); drawnow;

figure;
plot(Periods(1:t_plotmax), p_t(1:t_plotmax),'LineWidth',2,'Color',color1_ss); 
title('Price Path'); drawnow; 
orient landscape;

figure;
subplot(2,1,1);
plot(Periods(1:t_plotmax), rs_t(1:t_plotmax)*100,'LineWidth',2,'Color',color1_ss); hold on;
plot(Periods(1:t_plotmax), rb_t(1:t_plotmax)*100,'LineWidth',2,'Color',color1_shock); grid on;
title('Credit Market Rates'); ylabel('% Annualized'); legend('r^b','r^d');
grid on;
subplot(2,1,2);
plot(Periods(1:t_plotmax), rr_t(1:t_plotmax)*100,'LineWidth',2,'Color','k','LineStyle','-.'); hold on;
plot(Periods(1:t_plotmax), (rsp_t(1:t_plotmax)+rr_t(1:t_plotmax))*100,'LineWidth',2,'Color',color2_shock);
title('Policy Rates'); drawnow; grid on; ylabel('% Annualized'); legend('r^o^r','r^d^w')
% annotation('textarrow',[0.1,0.2],[0.9,0.6],...
%            'String','Ex-ante');
% annotation('textarrow',[0.1,0.2],[0.9,0.6],...
%            'String','Ex-ante');
if printit==1
    orient landscape;
    saveas(gcf,'F_mp_ratepath_temp','pdf');
end

figure;
plot(Periods, Y_t,'b'); hold on
plot(Periods, C_t,'r--');
plot(Periods, G_t,'k'); 
legend('Output','Consumption','Fed Profits'); 

figure;
subplot(2,1,1);
plot(Periods(1:t_plotmax), (Y_t(1:t_plotmax)/Y_ss-1)*100,'LineWidth',2,'Color',color1_ss); hold on;
title('Output - Deviation from Steady State'); ylabel('% from SS'); grid on;
axis tight; grid on;
subplot(2,1,2);
plot(Periods(1:t_plotmax), rr_t(1:t_plotmax)*100,'LineWidth',2,'Color','k','LineStyle','-.'); hold on;
plot(Periods(1:t_plotmax), (rsp_t(1:t_plotmax)+rr_t(1:t_plotmax))*100,'LineWidth',2,'Color',color2_shock);
title('Policy Rates'); drawnow; grid on; ylabel('% Annualized'); legend('r^o^r','r^d^w')
if printit==1
    orient landscape;
    saveas(gcf,'F_mp_ratepath_y_temp','pdf');
end

figure;
% plot(Periods, D_t,'b'); hold on
subplot(2,1,1);
plot(Periods(1:t_plotmax), (B_t(1:t_plotmax)/B_ss-1)*100,'LineWidth',2,'Color',color1_ss); hold on;
title('Borrowing and Lending Volume'); ylabel('% from SS'); grid on;
axis tight; grid on;
subplot(2,1,2);
plot(Periods(1:t_plotmax), rr_t(1:t_plotmax)*100,'LineWidth',2,'Color','k','LineStyle','-.'); hold on;
plot(Periods(1:t_plotmax), (rsp_t(1:t_plotmax)+rr_t(1:t_plotmax))*100,'LineWidth',2,'Color',color2_shock);
title('Policy Rates'); drawnow; grid on; ylabel('% Annualized'); legend('r^o^r','r^d^w')
if printit==1
    orient landscape;
    saveas(gcf,'F_mp_credit_temp','pdf');
end

figure;
% plot(Periods, D_t,'b'); hold on
subplot(2,1,1);
plot(Periods(1:t_plotmax), (p_t(1:t_plotmax)/p_ss-1)*100,'LineWidth',2,'Color',color1_ss); hold on;
title('Price Path'); ylabel('% from SS'); grid on;
axis tight; grid on;
subplot(2,1,2);
plot(Periods(1:t_plotmax), rr_t(1:t_plotmax)*100,'LineWidth',2,'Color','k','LineStyle','-.'); hold on;
plot(Periods(1:t_plotmax), (rsp_t(1:t_plotmax)+rr_t(1:t_plotmax))*100,'LineWidth',2,'Color',color2_shock);
title('Policy Rates'); drawnow; grid on; ylabel('% Annualized'); legend('r^o^r','r^d^w')
if printit==1
    orient landscape;
    saveas(gcf,'F_mp_price_temp','pdf');
end

figure;
plot(Periods(1:t_plotmax), (p_t(1:t_plotmax)/p_t(1)-1)*100,'LineWidth',2,'Color',color2_shock); title('Price Path'); drawnow;
grid on;

figure;
subplot(2,1,1);
plot(Periods, B_t,'LineWidth',2,'Color',color1_ss); hold on; grid on;

figure;
plot(s_vec, f_ss,'b.-'); hold on
plot(s_vec, f_0,'r--');
plot(s_vec, f_t(:,end),'k*');
legend('f_ss','f_0','f_T'); grid on;

figure
plot(Periods(1:t_plotmax),f_t(1,1:t_plotmax));

figure
surf(Periods(1:t_plotmax),s_vec(1:s_plotmax),f_t(1:s_plotmax,1:t_plotmax));
