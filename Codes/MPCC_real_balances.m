%% MPCC_implicit_general

% Main Code to Solve and Simulate the MPPC model
% (c) Bigio and Sannikov
%
% October 2016
%
% Relative to previous versions, now solves a general policy
% that is, it allows for discount window rates, rate on reserves,
% transfers/fed equity policy and money supply rule.


% Main block happens here
clear; 
close all;

%% Notes
% - Solves the model in real terms
% - No default.
% - Includes aggregate demand externality, transfers

%% Running Preferences
plotiter=0; % set to 1, if want to see algorithm converge
plotit  =0; % 
breakit =0; % 
printit =0; % 

%% Globals
MPCC_globals_MP; % Define list of globals

%% Code Parameters
% Approximate Amount of grid points
N = 160          ; % Number of Gridpoints in Space

% Time Parameters
T     = 80       ; % Vector of Intervals
Titer = 1        ; % Iterations

% From older Code
T_pre = 45       ; % Time Lapse prior to transition
T_post= 55       ; % Time Lapse after transition
Tspan = [1/12 0] ; % Vector Time, within interval

% State-Space Grid
dt = 1         ; % Time Dimention
dt_f = dt/Titer; % A step size for the KFE

% Tolerance 
tol=10e-5;
cond=2*tol;

options = optimset('TolFun',1e-3);

% Plot Properties
s_plotmax=round(N/2);
t_plotmax=80 ;

color1_ss   =[0.2 0.2 0.6];
color1_shock=[0.8 0.2 0.2];
color2_shock=[0.6 0.2 0.2];

%% Model Parameter Definitions
% [-] Pass Labels to Functions

% Auxiliaries
iter=0;

% Model Parameters - At Steady State
rho = 0.05  ; % agent's individual discount rates -> presumably need rho >= rs to make sure savings don't blow up
gamma = 2   ; % agent's risk aversion
w1 = 2      ; % mean low process
w2 = 10     ; % mean high process 
s1 = 0      ; % volatility low process
s2 = 10     ; % volatility 

% Borrowing Limit - nominal
s_bar  = -0*w2   ;  % debt limit
s_bl   = s_bar   ;  % borrowing limit - below this point, you can accumulate interest but not borrow

% Vectorizations
mu_w    = [w1 w2]';
sigma_w = [s1 s2]';

% Fee to pay back
fee = 2*abs(s_bar);

% Bliss point
c_bliss=-w1*0; % Potentiall bliss point to break homogeneity assumption - useful to undo risk-aversion effect
% -> code not ready, we need to modify this assumption

% Real transers
tau=0.00*w1;

%% Functional Forms 
if gamma~=1
    URF=@(c) ((c-c_bliss).^(1-gamma)-1)/(1-gamma); % Utility Return Function
    Upinv=@(c) (c-c_bliss).^(-1/gamma); % Utility Return Function
else
    URF=@(c)   log(c-c_bliss); % Utility Return Function
    Upinv=@(c) exp(c-c_bliss); % Utility Return Function
end

%% Policy and Assumed Prices
rd = 0.00   ; % Deposit Rate
rb = 0.00   ; % Borrowing Rate
rs = 0.00   ; % Savings Rate
p  = 1      ; % Initial Steady-State Price Level 

%% Deduced Parameters 
mpc    = (rd + (rho - rd)/gamma)        ; 
A      = 1/(rho-(1-gamma)*(rd-mpc))     ; 
s_max = 1.2*w2/rho                      ;  % consider savings up to 4 times PV of expected risky income 

% Borrowing Limit - Static
c_dl = rs*s_bar/p + w1; % Stays Constant

% Borrowing Limit - Dynamic
c_bl(1)  = c_dl  ; 
c_bl(2:N,1)= w2  ; 

% Idea is that below s_bar_cash, you can accumulate debt but not borrow.
% this means that debt can only increase at r_rate, and not faster rate
% mu_s>rs so (w2 > c(2:end)) in that region

%% Discretization of State Space
% Savings space and consumption space Gus
s_vec  = linspace(s_bar, s_max,N)'           ; % Space of (nominal) savings 

% Initialized Variables
mu_ss  = zeros(1,N); % Drift at Steady State
sigma_ss = zeros(1,N); % Volatility at Steady State

%% Steady State Values
rs_ss=rs;

%% Experiments - Determinisitic paths
% Only allowing Technology / policies to change 

% [1.a ] *Means and Volatilities along Transition
shock_mu   =1                               ;
mu_w_t     =mu_w*ones(1,T)                  ;
mu_w_t(2,T_pre+1:T_post)= mu_w(2)*shock_mu  ;

% [1.b ]Means and *Volatilities along Transition
sigma_w_t=sigma_w*ones(1,T)  ;

% path
shock_sigma=1.0;
sigma_w_t(2,T_pre+1:T_post)= s2*shock_sigma;

% [2] Rate in Vector of Time and space
r_t(1:T)                   = rs_ss ;
rd_t(1:T)                  = rd;
rs_t(1:T)                  = rs;
rb_t(1:T)                  = rb;

% Add Credit Variables
rb_vec(1:T_pre)              = rb;
rb_vec(T_pre+1:T_pre+T_post) = rb;
rb_vec(T_pre+T_post+1:T)     = rb;

% [3] Borrowing Limits 
shock_s_bl=0.0;
s_bl_t(1:T)=s_bl;
s_bl_t(T_pre+1:T_post)=s_bl*shock_s_bl;

% [4] Transfer Shocks - Expected
shock_T=0.0; % Rate of money growth do to transfers 
T_t(1:T)=0  ;
T_t(T_pre+1:T_post)=shock_T*dt_f;

% [5] Helicopter Drop
Delta_ingrid=5                      ; % Unexpected shock

% write down paths
paths.r_t=r_t            ;
paths.rs_t=rs_t          ; 
paths.rb_t=rb_t          ; 
paths.rd_t=rd_t          ; 
paths.mu_w_t=mu_w_t      ;
paths.sigma_w_t=sigma_w_t;
paths.s_bl_t=s_bl_t      ;
paths.T_t=T_t            ; 

% Paths for exogenous 
paths_list={'r_t','rs_t','rb_t','rd_t','mu_w_t',...
    'sigma_w_t','s_bl_t','T_t'};
paths_list_tit={'r_t','rs_t','rb_t','rd_t','\mu_2',...
    '\sigma_2','s^b^l__t','T_t'};

%% Add Separate code to describe shocks
MPPC_report_shocks

%% [1] Solving Initial Steady State
% Initial conjecture for consumption
c_vec = c_bliss+linspace(w1, w2, N)'+rs_ss*s_vec      ; % Initial guess - permanent income hypothesis

%Solution of HJB Equation using Policy Iteration Method
r=rs_ss;
c=c_vec;
% s=s_vec;

% Update credit limit to steady-state value
c_bl(1)  = c_dl       ; % 
c_bl(2:N,1)= w1       ; % Forces a switch to safe technlogy alternatively

% timers
tic
iter=0;

% Solve for steady state
MPCC_solve_real_steady;

% Solving for Steady-State Distributions
[f_ss]=KFE_ss_implicit(mu_ss,sigma_ss,s_vec,dt);
f_ss=f_ss'                  ; 

% Real Aggregate Endowment and Consumption 
Y_ss=sum(f_ss(:)*w2) + f_ss(1)*(w1 - w2);
C_ss=f_ss'*c_ss                         ;
GoodsZ=(Y_ss-C_ss)/Y_ss*100                 ;
RS_ss=f_ss'*s_vec                       ; % Real Savings

% Register Initial Conditions
steady.V_ss1     = V_ss       ; 
steady.mu_ss1    = mu_ss      ;
steady.sigma_ss1 = sigma_ss   ;
steady.f_ss1     = f_ss       ;
steady.C_ss1     = C_ss       ;
steady.Y_ss1     = Y_ss       ;
steady.RS_ss1    = RS_ss      ;
% steady.pi_ss1    = pi_ss      ;

% Steady State Inflation - 
% Government Budget Constraint: 
% dM_t= p_t*tau*d_t;
% M_t = p_t*RS_ss   ; %-> quantity equation yields:
% dM_t/M_t= RS_ss^(-1)*tau*d_t;
% but dM_t = dp_t*RS_ss -> dp_t/p_t=RS_ss^(-1)*tau*dt
pi_ss=tau/RS_ss              ;

% Money Stock
p_0 = 1                        ;
M0_0=  RS_ss                   ;
p_t = p_0*exp(pi_ss*dt*(1:T))  ;

% Nominal Variables - Here....
% MStock_ss1=MStock_ss         ;
% steady.MStock_ss1=MStock_ss  ;
% Registering Values
% MStock_0=MStock_ss1*(1+Delta_M)    ;

%% Unexpected time zero shock
% Equivalent to translation to the right
% One time Transfer - Translation Shock
trans=(s_vec(Delta_ingrid)-s_vec(1)); % Compute real transfer size
Delta_R=(trans)/RS_ss               ; % Real Transfer

% Cummulative Distribution
F_ss = cumsum(f_ss);

% New Distribution afte unexpected shock:
f_0=zeros(N,1)                                      ;
f_0(Delta_ingrid)=f_ss(1)                           ;
F_0=zeros(N,1)                                      ;
F_0(Delta_ingrid)=f_ss(1)                           ;
F_0(Delta_ingrid+1:end)= interp1(s_vec,F_ss,(1+r-Delta_R)^(-1)*s_vec(Delta_ingrid+1:end)-trans,'linear');
F_0(isnan(F_0))=1;
f_0(Delta_ingrid+1:end)= F_0(Delta_ingrid+1:end)-F_0(Delta_ingrid:end-1);
 
% Initial Condition
init_con.f_0=f_0;
f_pre=f_ss;

% First Conditino 
Y_0=sum(f_0(:)*w2) + f_0(1)*(w1 - w2);
C_0=f_0'*c_ss                         ;
Z_Y_0=(Y_ss-C_ss)/Y_ss*100                 ;
RS_ss=f_0'*s_vec                       ; % Real Savings


%% Define final steady state
% p_ss=p;
% 
% % Register Initial Conditions
steady.V_ss     = V_ss       ; 
steady.mu_ss    = mu_ss      ;
steady.sigma_ss = sigma_ss   ;
steady.f_ss     = f_ss       ;
steady.C_ss     = C_ss       ;
steady.Y_ss     = Y_ss       ;
steady.RS_ss    = RS_ss      ;
steady.pi_ss   = pi_ss      ;


%% Plots for Comparisons among Methos.
% Pictures
if plotit==1
    figure;
    hold on;
    
    % Initial Value 
    plot(s_vec(2:end),f_ss(2:end),'Color',color1_ss,'LineWidth',2);
    % h2=area(s_vec,f_ss,'LineWidth',1,'FaceColor',color2_shock,'FaceAlpha',0.6); 
    line([s_vec(1) s_vec(1)],[0 f_ss(1)],'Color',color1_ss,'LineWidth',2);
    scatter(s_vec(1),f_ss(1),'MarkerFaceColor',color1_ss);
    title('Invariant Distribution of Real Balances'); axis tight; xlabel('Real Balances'); ylabel('Density');
    grid on;

    % leyend
%     leyenda=[h2(1)];
   % legend(leyenda,'Initial Distribution');

    if printit==1
        orient landscape
        saveas(gcf,'F_rb_f_ss','pdf');
    end
    
    % Final Figure
    plot(s_vec(Delta_ingrid+1:end),f_0(Delta_ingrid+1:end),'Color',color1_shock,'LineWidth',2);
   %  h1=area(s_vec(Delta_ingrid:s_plotmax),[f_0(Delta_ingrid+1) f_0(Delta_ingrid+1:s_plotmax)'],'LineWidth',1,'FaceColor',color1_shock,'FaceAlpha',0.6); 
    line([s_vec(Delta_ingrid) s_vec(Delta_ingrid)],[0 f_0(Delta_ingrid)],'Color',color1_shock,'LineWidth',2);
    scatter(s_vec(Delta_ingrid),f_0(Delta_ingrid),'MarkerFaceColor',color1_shock);
    title('Invariant Distribution of Real Balances'); axis tight; xlabel('Real Balances'); ylabel('Density');
    grid on;
        
%     x = [s_vec(1) s_vec(Delta_ingrid)]      ;
%     y = [f_0(Delta_ingrid) f_0(Delta_ingrid)]/2/f_0(Delta_ingrid)   ;
%     [figx figy] = dsxy2figxy(x,y)                 ;
%      annotation('textarrow',figx ,[1/3 1/3],'String','Injection');
% 
%     x = [1/6 1/3]*s_vec(end)      ;
%     y = [f_0(Delta_ingrid) f_0(Delta_ingrid)]/2/f_0(Delta_ingrid)   ;
%     [figx figy] = dsxy2figxy(x,y)                 ;
%     annotation('textarrow',figx ,[1/3 1/3],'String','highher mass');
%      
%     leyenda=[h1(1) h2(1)];
%     legend(leyenda,'After Shock','Initial Distribution');


%     figure
%     plot(s_vec,steady.V_ss1 ,'b'); hold on;
%     plot(s_vec,steady.V_ss2b,'r--');
%     plot(s_vec,steady.V_ss2,'g-.');
%     grid on; legend('before','after'); title('Value of Money Before and After HD');
%     
%     figure
%     plot(s_vec,steady.mu_ss1 ,'b'); hold on;
%     plot(s_vec,mu_ss,'r--'); hold off;
%     grid on; legend('before','after'); title('Value of Money Before and After HD');
%     
%     figure
%     plot(s_vec,steady.sigma_ss1 ,'b'); hold on;
%     plot(s_vec,sigma_ss,'r--'); hold off;
%     grid on; legend('before','after'); title('Value of Money Before and After HD');
%     
%     figure
%     h3=area(s_vec(1:s_plotmax),[f_ss2(2) f_ss2(2:s_plotmax)'],'LineWidth',1,'FaceColor',color1_ss,'FaceAlpha',0.9); hold on;
%     line([s_vec(1) s_vec(1)],[0 f_ss2(1)],'Color',color1_ss,'LineWidth',2);
%     scatter(s_vec(1),f_ss2(1),'MarkerFaceColor',color1_ss);
% 
%     h2=area(s_vec(1:s_plotmax),[steady.f_ss2a(2) steady.f_ss2a(2:s_plotmax)'],'LineWidth',1,'FaceColor',color1_shock,'FaceAlpha',0.9); 
%     line([s_vec(1) s_vec(1)],[0 steady.f_ss2a(1)],'Color',color1_shock,'LineWidth',2);
%     scatter(s_vec(1),steady.f_ss2a(1),'MarkerFaceColor',color1_shock);
% 
%     h1=area(s_vec(1:s_plotmax),[steady.f_ss2(2) steady.f_ss_2b(2:s_plotmax)'],'LineWidth',1,'FaceColor',color2_shock,'FaceAlpha',0.9); 
%     line([s_vec(1) s_vec(1)],[0 steady.f_ss2(1)],'Color',color2_shock,'LineWidth',2);
%     scatter(s_vec(1),steady.f_ss2(1),'MarkerFaceColor',color2_shock);
%     
%     title('Invariant Distribution of Nominal Wealth'); axis tight;
%     grid on;
    
end

%% Guess and Solve p_t
% Guess for p_t;
% p_t_o = p_ss*ones(1, T-1)    ; % Price in each Period 
pi_t_o = [linspace(-0.09*Delta_R,0,20) zeros(1,60)];

% [*] Solving for the price path
res_sup=@(pi_t) MPCC_path_solve_rb(pi_t,init_con,steady,paths);
tic
pi_t_sol=fsolve(res_sup,pi_t_o,options);
disp('Time to Find Solution');
toc
resY_out=res_sup(pi_t_sol);

%% Recomputing Equilibrium with Price Solution
% p_t=[p_t_sol steady.p_ss];
pi_t=[pi_t_sol];
    
% toc
% pi_t = linspace(-0.09*Delta_R,0,T);
MPCC_path_rb

% Nominal Values
p_t=exp(cumsum(pi_t));

%% Evolution of Real Wealth
% rf_t=f_t./(ones(N,1)*p_t);

% Compute Solution
Periods = (1:1:T);
if plotit==0  
    figure;
    plot(Periods, Z_y,'b'); title('Excess Supply of Goods'); drawnow;

    % Plot Price Path
    figure;
    extra_t=(-10:1:0);
    plot([0 Periods],[NaN p_t],'b-'); title('Price Path'); hold on;
    scatter(0,1/(1+Delta_R),'b','filled'); title('Price Path'); drawnow;
    plot(extra_t,extra_t*0+p/(1+Delta_M),'b-.'); title('Price Path'); hold on;
    scatter(1,p_t(1),'b'); title('Price Path'); 
    scatter(Periods(end),p_ss,'r','filled'); title('Price Path'); drawnow;
    legend('pre-shock path','initial price','after-shock path','after-shock jump','new steady state'); grid on; xlim([extra_t(1) Periods(T)]);
    if printit==1
        orient landscape
        saveas(gcf,'F_rb_p_t','pdf');
    end

    % Plot Output Path
    y_ss=Y_t(end);
    y_dev_t=(Y_t./y_ss(end)-1).*100;
    y_ss_dev=(Y_ss/y_ss(end)-1)*100;
    figure;
    extra_t=(-10:1:0);
    plot([0 Periods],[NaN y_dev_t],'b'); title('Output Path'); hold on;
    scatter(0,0,'b','filled'); 
    plot(extra_t,extra_t*0+0,'b'); 
    scatter(1,y_dev_t(1),'b'); drawnow;
    scatter(Periods(end),y_ss_dev,'b','filled'); drawnow;
    legend('path','initial output'); grid on; xlim([extra_t(1) Periods(T)]);
    if printit==1
        orient landscape
        saveas(gcf,'F_rb_y_t','pdf');
    end

    % Plot Wealth Distribution
    figure
    surf(Periods(1:t_plotmax),s_vec(1:s_plotmax),f_t(1:s_plotmax,1:t_plotmax)); axis tight;
    xlabel('time'); ylabel('Nominal Balances'); title('Nominal Wealth Distribution Evolution');
    if printit==1
        orient landscape
        saveas(gcf,'F_rb_f_t','pdf');
    end
end

%% Welfare

%% Debt-Deflation

% figure(pic+5)
% s_plotmax=120;
% surf(Periods,s_vec(1:s_plotmax),rf_t(1:s_plotmax,:)); axis tight;
% xlabel('time'); ylabel('chips'); title('Real Wealth Distribution Evolution');