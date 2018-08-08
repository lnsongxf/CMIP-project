%% MPCC_implicit_general

% Main Code to Solve and Simulate the MPPC model
% (c) Bigio and Sannikov
%
% October 2016
% Main block happens here
clear all; 
close all;
nameplot = 'rp_ext2';

addpath(genpath([cd '\Codes']));
% Graphs
path_g       = [cd '\Graphs'];
imprime      = @(x) print( gcf, '-depsc2', [path_g filesep x]);
imprpdf      = @(x) eps2pdf( [path_g filesep x '.eps']);
formataxis   = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 19, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
formatlegend = @(x) set(legend, 'Location', x, 'Orientation', 'Vertical', 'Box', 'On', 'Fontsize', 20, 'Fontangle', 'normal');
label_x      = @(x) xlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15,'interpreter','latex');
label_y      = @(x) ylabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15,'interpreter','latex');
label_z      = @(x) zlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15,'interpreter','latex');
otitle       = @(x) title(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15,'interpreter','latex');

%% Definition of Global Variables
MPCC_globals; % Define list of globals

flag_unant =0;  % [1] Unanticipated shock
cc         =1;  %figure counter
solveit    =1;  % solveit:
                % [0] uses bisection algorithm to find interest rate
                %     path.
                % [1] uses LM algorithm and solver function
                

%% [V] Policy and Assumed Prices
mpregime='RP' ;

% set to:
% -------------------------------------------------------------------------
%
% FP: Fisherian principle
% RP: Ricardian principle
% BP: Bagehot principle
% -------------------------------------------------------------------------               
                
%% Running Preferences
plotiter  =0; % set to 1, if want to see algorithm converge
plotit    =0; % 
breakit   =0; % 
printit   =1; % set to 1, if you want to save figures.
plotib    =0; % Interbank market plots
flag_compute = 1;

%% [I] Code Parameters and shocks parameters
% Convergence Criterion
tol_dyn = 1e-3;
xi      = 1e-5;
tol     = 1e-6;
cond    = 2*tol;
options = optimset('TolFun',1e-8,'Display','iter'); % At Optimization

% Approximate Amount of grid points
N     = 500;            % Number of Gridpoints in Real Wealth Space

% Time Parameters
T     = 100      ;  % Time Horizon
Titer = 1        ;  
Delta = 10       ;  % Delta for solving Bellman equation at the stationary equilibrium.

% Periods for Transition
if flag_unant
    T_pre=1;
    T_post=11;
else
    T_pre = 15       ; % Time Lapse prior to transition
    T_post= 30       ; % Time Lapse after transition
end

% State-Space Grid
dt     = 1;

T      = T/dt;
T_pre  = T_pre/dt;
T_post = T_post/dt;
dt_f   = dt/Titer; % A step size for the KFE

% Shocks
shock_mu   = 1.0                            ;
shock_sigma= 1.0                            ;
rsp_shock  = 1                              ;
shock_s_bl = 0.8                            ;
shock_T    = 0                              ;  

% Clearing Condition
clearcond='Y'; % Y for goods market and S for asset market

%% [II] Plot Preferences
% Plot Properties
s_plotmax=round(N/4);
t_plotmax= 80;

% Steady State - Shock Colors
color1_ss   =[0.2 0.2 0.6];
color1_shock=[0.8 0.2 0.2];
color2_shock=[0.6 0.2 0.2];

%% [III] Model Parameter Definitions

% Model Parameters - Preferences and Technology
gamma = 0.8   ; % agent's risk aversion
rho   = 0.04 ; % agent's individual discount rates 
w1    = 50    ; % mean return - low intensity technology
w2    = 100   ; % mean return - high intensity technology
s1    = 0     ; % volatility - low intensity technology
s2    = 100   ; % volatility - high intensity technology
eta_l = 3/4   ; % 1-eta_l measure externality
N1    = 20    ;
N2    = 60    ;

% Initial Guess for interest rate rate
rs_o   = 0.03;

% Borrowing Limit - nominal
s_bar      = -10*w2        ;  % debt limit 
s_bl       = s_bar        ;  % index where 

% Vectorizations
mu_w    = [w1 w2]';
sigma_w = [s1 s2]';

% Bliss point
c_bliss = 0; 

% FED Portfolio block
Ef_ss       =     -0.2 ; % Governemnt's net worth/Deposits   
delta_ef    =     0.9  ; % adjustment of FED policy
delta_trans =     0.9  ; % speed of transfer adjustment

% Interbank Market block 
varrho  = 0.5 ; % Reserve Requirement
barlam  = 1.5 ; % Efficiency interbank market
eta     = 0.5 ; % Bargaining
omega   = 0.3 ; % Average size of shock - interbank

% Steady State Targets:
imonmu_ss = varrho;
i_m_ss    = 0.00;
isp_ss    = 0.02;   

% Target for real spread
rsp_ss     = 0.00;

% For Plots
step       = 0.0001;
surplus_cut=(varrho-omega*(1-varrho));
deficit_cut=(varrho+omega*(1-varrho));
delta_cut=deficit_cut-surplus_cut;
mu_vec=[(surplus_cut-0.005*delta_cut:step:surplus_cut)...
    (surplus_cut+step:step:deficit_cut-step)...
    (deficit_cut:step:deficit_cut+0.01*delta_cut)];

M_aux     = 100        ;
fedeqaux  = 000        ;

poldaux = isp_ss; % typical policy spread
r_m_aux = i_m_ss;   % this will be adjusted 

%% [IV] Functional Forms 
% Utility Functions
if gamma~=1
    URF=@(c) ((c-c_bliss).^(1-gamma)-1)/(1-gamma); % Utility Return Function
else
    URF=@(c)   log(c-c_bliss); % Utility Return Function
end

Upinv=@(dv) dv.^(-1/gamma)+c_bliss;              % Utility Return Function

%% [VI] Discretization of State Space
s_max = (w2/rho)*3.2                      ; 

% Savings space and consumption space Gus
s_vec  = linspace(s_bar, s_max,N)'        ; 
s_bl_index = (s_vec<=s_bl);

% Compute Deposits
index_o=find(s_vec>=0,1,'first')          ;

% Initialized Variables
mu_ss    = zeros(1,N); % Drift at Steady State
sigma_ss = zeros(1,N); % Volatility at Steady State

%% [IX] Solving Initial Steady State
% Set timers
tic
iter=0;

% Solve using unit price
p=1       ;
p_ss=1    ;
inputs.p=p; % Normalization of price level

% Initial conjecture for consumption
V = s_vec;

%Solution of HJB Equation using Policy Iteration Method
% c=c_vec;
inputs.V=V;

plot_zA=0;
if plot_zA==1   % upper bound lower than rho
        fplot(@(x) MPCC_solve_steady(x,inputs,paths),[-0.01 rho-1e-2]); grid on;
end

% Define and solve steady state rates
eqsteady=@(rs_o) MPCC_solve_steady(rs_o,inputs);
rs_ss=fsolve(eqsteady,rs_o,options);

% Update Rate 
rs=rs_ss                    ;

%-----------------------------------
% Recomputing solution at steady state
%-----------------------------------
% % Construct Interest Rate Vector
% [1] Solve for Steady State for each case
% Construct Interest Rate Vector
r_vec=s_vec       ; % Initialize it
r_vec(s_vec>0)=rs ;
if any(strcmp(mpregime,{'FP', 'RP'}));
    r_vec(s_vec<=0)=rs; % pin down value
elseif strcmp(mpregime,{'BP'})
    r_vec(s_vec<=0)=rs+rsp_ss; 
end

% Compute the Steady State
MPCC_steady;

disp('*****************************************');
disp('** Achieved Convergence with Clearing****');
disp('*****************************************');
disp(['Clearing in Goods: ' num2str(Z_Y_ss)]);
disp(['Clearing in Assets: ' num2str(Z_S_ss)]);
disp('*****************************************');
disp('*********** Beggining Transition** ******');
disp('*****************************************');

% Value of Real Transfers
TE_ss         = Ef_ss*D_ss          ;
RT_ss         = TE_ss*rs            ;
FT_ss         = rsp_ss*B_ss         ;

% mp_implementation;

% Saving steady state variables and initial conditions
% Recording Terminal Steady State
steady.V_ss      = V_ss      ; 
steady.sigma_ss  = sigma_ss  ;
steady.f_ss      = f_ss      ;
steady.C_ss      = C_ss      ;
steady.Y_ss      = Y_ss      ;
steady.c_ss      = c_ss      ;
steady.muF_ss    = muF_ss    ;
steady.muB_ss    = muB_ss    ;
steady.D_ss      = D_ss      ;
steady.B_ss      = B_ss      ;
steady.TE_ss     = TE_ss     ; 
steady.RT_ss     = RT_ss     ; % Value of real transfers
steady.FT_ss     = FT_ss     ; % Value of real transfers
steady.rs_ss     = rs_ss     ; 
steady.mu_ss     = mu_ss     ;
steady.y1_ss     = y1        ;
steady.y2_ss     = y2        ;
steady.rsvec=r_vec.*s_vec;

% Initial conditions - distribution
f_0          = f_ss;
init_con.f_0 = f_0;

%% [VII] Experiments - Determinisitic paths

% [1.a ] *Means and Volatilities along Transition
mu_w_t     = mu_w*ones(1,T)                 ;
mu_w_t(2,T_pre+1:T_post)= mu_w(2)*shock_mu  ;

% [1.b ]Means and *Volatilities along Transition
sigma_w_t=sigma_w*ones(1,T)                ;
sigma_w_t(2,T_pre+1:T_post)= s2*shock_sigma;

% [2.a] Vector of Policy Rates - interest on reserves
rr_shock                     = i_m_ss   ;
rr_t(1:T)                    = i_m_ss   ; 
rr_t(T_pre+1:T_post)         = rr_shock ;

% [2.a] Vector of Policy Rates - interest on reserves
rsp_t(1:T)                    = rsp_ss    ; % Interest on reserves
rsp_t(T_pre+1:T_post)         = rsp_shock*rsp_ss;
% rsp_t(T_post+1:T_post+5)      = [0.6 .7 .8 .9 .95]*rsp_ss;

% [3] Borrowing Limits 
temp       = s_bl*shock_s_bl          ;
[~,inds]   = min(abs(bsxfun(@minus,s_vec,temp)))   ;
s_bl_t(1:T)=s_bl                      ;
s_bl_t(T_pre+1:T_post)= s_vec(inds)   ;
s_bl_index_t=zeros(T,1)               ;
for tt=1:T
    s_bl_index_t(tt)= sum(s_vec<=s_bl_t(tt))    ;
end

% [4] Transfer Shocks - Expected
T_t(1:T) = 0                       ;
T_t(T_pre+1)=shock_T               ; 
if T_post-T_pre>2
    for tt=T_pre+2:T_post
        T_t(tt)=(T_post+1-tt)/(T_post-T_pre+1)*shock_T ;
    end
end

% [5] Helicopter Drop
Delta_ingrid=1                      ; % Unexpected shock

% write down paths
paths.rr_t      = rr_t       ;
paths.rsp_t     = rsp_t      ;
paths.mu_w_t    = mu_w_t     ;
paths.sigma_w_t = sigma_w_t  ;
paths.s_bl_t    = s_bl_t     ;
paths.T_t       = T_t        ; 
paths.s_bl_index_t=s_bl_index_t;

% Paths for exogenous 
paths_list={'rr_t','rsp_t','mu_w_t','sigma_w_t','s_bl_t','T_t'};
paths_list_tit={'r^o^r_t','r^d^w_t','\mu_2','\sigma_2','s^b^l_t','T_t'};

%% [VIII] Run a code that reports the shocks
MPPC_report_shocks;

%% [X] Solving Transition 
% Guess for prices
p_t    = ones(T,1)                     ;
rs_t_o = rs_ss*ones(1, T)              ;

% Plot Fed Experiment for Guess
T_t=paths.T_t                          ;
index_T=(T_t~=0)                       ;
index_T_T=find(index_T==1,1,'last')    ;

% Solving for the price path
switch mpregime
    case {'FP'}
        eqpath=@(x)      MPCC_path_solve_FP(x,init_con,steady,paths)     ;
    case {'RP'}
        eqpath=@(rs_t_o) MPCC_path_solve_RP(rs_t_o,init_con,steady,paths);
    case {'BP'}
        eqpath=@(rs_t_o) MPCC_path_solve_BP(rs_t_o,init_con,steady,paths);
end

% Solving the main function
if flag_compute
    tic
    if solveit==1
        options_dyn=optimset('Algorithm',{'levenberg-marquardt',0.01},'MaxFunEvals',12000,'MaxIter',400,'TolFun',1e-5,'Display','iter');
        rs_t=fsolve(eqpath,rs_t_o,options_dyn);
    else
        SOLVER_MA;
        rs_t = r_t;
    end
    disp('Time to Find Solution');
    toc
else 
    load(['DD' nameplot '.mat']);
    rs_t = endopath.rs_t;
end
    
% Steady State Values
V=steady.V_ss;
if ~isempty(Dr_zlb)
    rsp_t=max(rsp_t,Dr_zlb);
end
switch mpregime
    case {'FP'}
        rb_t=rs_t;
    case {'RP'}
        rb_t=rs_t;
    case {'BP'}
        rb_t=rs_t+rsp_t;
end

% Backout Transfer set
T_t=paths.T_t                          ;
index_T=(T_t~=0)                       ;
index_T_T=find(index_T==1,1,'last')    ;
RT_t=RT_ss*ones(1,T)                   ;
RT_t_aux=RT_t                          ;
TE_t=TE_ss*ones(1,T)                   ;
TE_t_aux=TE_t                          ;

% Equity Path
for tt=1:T
    if tt<=T_post
        RT_t(tt  ) = (1+T_t(tt))*RT_ss                  ;
        TE_t(tt+1) = (1+rs_t(tt))*TE_t(tt)-RT_t(tt)   ;
    elseif tt>T_post
        aux_weight=delta_trans^(tt-T_post)                    ;
        TE_t_aux(tt+1) = TE_ss+delta_ef*(TE_t(tt)-TE_ss)      ;
        RT_t_aux(tt) = (1+rs_t(tt))*TE_t(tt)-TE_t_aux(tt+1)   ; % Find transfer rule that clears things
        RT_t(tt)     = (1-aux_weight)*RT_t_aux(tt)...
            +aux_weight*RT_t(T_post);
        TE_t(tt+1) = (1+rs_t(tt))*TE_t(tt)-RT_t(tt) ;
    end
end

% Solve again - with the new interest rate
MPCC_path                   ;

% mass of constraint agents
fc_t = zeros(T,1);  
for t=1:T
    temp    = paths.s_bl_index_t(t); 
    fc_t(t) = sum(f_t(1:temp,t)*ds);
end

% Save
endopath.rs_t=rs_t;
endopath.rb_t=rb_t;
endopath.TE_t=TE_t;
endopath.RT_t=RT_t;
endopath.Y_t=Y_t;
endopath.B_t=B_t;
endopath.D_t=D_t;
endopath.f_t=f_t;

eval(['save DD' nameplot ' c_t V_t f_t D_t B_t s_vec rs rs_t rs_ss f_ss steady paths endopath']);
rs_t(1) = rs_ss;
rb_t(1) = rs_ss+rsp_ss;
MPCC_nominalimplementation;

%% [XI] Welfare Calculations
% Ex-Ante Welfare
V_ea=V_t(:,1);
if gamma~=1
    CE_ea=((1-gamma)*rho*V_ea+1).^(1/(1-gamma))+c_bliss;
else
    CE_ea = exp(rho*V_ea)+c_bliss;
end

EV_ea=f_t(:,1)'*V_ea*ds;
ECE_ea=f_t(:,1)'*CE_ea*ds;
ECE_ea_Y=(1-ECE_ea/steady.Y_ss)*100;

% Steady-State Value
if gamma~=1
    CE_ss=((1-gamma)*rho*steady.V_ss+1).^(1/(1-gamma))+c_bliss;
else
    CE_ss = exp(rho*steady.V_ss)+c_bliss;
end

E_V_ea=f_ss(:)'*steady.V_ss*ds;
ECE_ss=f_ss(:)'*CE_ss*ds;
ECE_ss_Y=(1-ECE_ss/steady.Y_ss)*100;

%% [XII] Main Plots
Periods = (1:1:T)*dt;
% Plot Values
figure(cc);
plot(CE_ea,f_t(:,1),'LineWidth',4); hold on;
plot(CE_ss,f_ss,'r--','LineWidth',4); 
label_x('Certainty Equivalent'); label_y('Distribution'); grid on;
legend('Ex-Ante','Steady-State');
formataxis(gca); axis tight;
h=text(CE_ss(10),f_ss(10),['\leftarrow CE % Output Loss: ' num2str(ECE_ea_Y-ECE_ss_Y,3)]);
set(h,'fontsize',18);
formatlegend('Best');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

% Figures for the Draft
figure(cc);
plot(Periods(1:t_plotmax)-1, rs_t(1:t_plotmax)*100,'LineWidth',4,'Color',color1_ss); hold on;
plot(Periods(1:t_plotmax)-1, rb_t(1:t_plotmax)*100,'--o','LineWidth',4,'Color',color1_shock); 
plot(Periods(1:t_plotmax)-1, rs_ss*100+0*rb_t(1:t_plotmax)*100,'k:','LineWidth',1); grid on;
axis tight; [yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
legend('r^a','r^l','r^{a}_{ss}');
label_x('time'); label_y('$\%$ annual');
% otitle('Real Borrowing/Lending Rates');
formatlegend('Best'); formataxis(gca); 
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

% Fed Policy
figure(cc)
plot(Periods(1:t_plotmax)-1, TE_t(Periods(1:t_plotmax))./Y_t(Periods(1:t_plotmax))*100,'LineWidth',4,'Color',color1_ss); hold on;
plot(Periods(1:t_plotmax)-1, TE_ss/Y_ss*100+0*TE_t(Periods(1:t_plotmax)),'k:','LineWidth',1); grid on;

[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
% otitle('$\mathcal{E}_{t}/Y_{t}$ deviation from $\mathcal{E}_{ss}/Y_{ss}$'); 
formataxis(gca); axis tight;
label_x('time'); label_y('$\%$');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc)
plot(Periods(1:t_plotmax)-1, RT_t(Periods(1:t_plotmax))./Y_t(Periods(1:t_plotmax))*100,'LineWidth',4,'Color',color1_ss); hold on;
plot(Periods(1:t_plotmax)-1, RT_ss/Y_ss*100+0*RT_t(Periods(1:t_plotmax)),'k:','LineWidth',1); grid on;
[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
% otitle('$T_{t}/Y_{t}$ deviation from $T_{ss}/Y_{ss}$'); 
formataxis(gca); axis tight;
label_x('time'); label_y('$\%$');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc);
plot(Periods(1:t_plotmax)-1, (Y_t(1:t_plotmax)/Y_ss-1)*100,'LineWidth',4,'Color',color1_ss); hold on;
plot(Periods(1:t_plotmax)-1, 0*rb_t(1:t_plotmax)*100,'k:','LineWidth',1); grid on;
axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
label_x('time'); label_y('$\%$');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc);
plot(Periods(1:t_plotmax)-1, (B_t(1:t_plotmax)/B_ss-1)*100,'LineWidth',4,'Color',color1_ss); hold on; 
plot(Periods(1:t_plotmax)-1, (D_t(1:t_plotmax)/D_ss-1)*100,'--o','LineWidth',4,'Color',color1_shock);  grid on;
axis tight; [yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
% otitle('$B_{t},A_{t}$ deviation from $B_{ss},A_{ss}$'); 
legend('B_t','A_t');
formataxis(gca);
formatlegend('Best');
label_x('time'); label_y('$\%$');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc);
plot(s_vec, f_ss,'b-'); hold on
plot(s_vec, f_t(:,T_pre+1/dt));
plot(s_vec, f_t(:,end),'k');
legend('f_ss','f_t(init-shock)','f_T'); grid on;

if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc)
plot(Periods(1:t_plotmax)-1,fc_t(1:t_plotmax)*100,'linewidth',4); grid on;
axis tight; [yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
% otitle('Mass at the constraint'); 
formataxis(gca); 
label_x('time'); label_y('$\%$'); 
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc)
surf(Periods(1:t_plotmax)-1,s_vec(2:s_plotmax),f_t(2:s_plotmax,1:t_plotmax),'edgealpha',.1);
hold on;
scatter3(Periods(1:t_plotmax)-1,s_vec(1)*ones(1,t_plotmax),f_t(1,1:t_plotmax)); 
axis tight;
% otitle('Distribution of Wealth');
formataxis(gca); set(gca,'box','off');
label_x('time'); label_y('real wealth');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

% Tests
figure(cc)
h=surface(Periods(1:t_plotmax)-1,s_vec(2:s_plotmax),f_t(2:s_plotmax,1:t_plotmax),'edgealpha',.01);
alpha(h,0.7);
% otitle('Distribution of Wealth'); 
label_x('time'); label_y('real wealth');
axis tight; 
formataxis(gca); set(gca,'box','off');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end    
cc=cc+1;

% Tests
figure(cc)
h=mesh(Periods(1:t_plotmax)-1,s_vec(2:s_plotmax),f_t(2:s_plotmax,1:t_plotmax));
alpha(h,0.5); 
axis tight; grid on;
% otitle('Distribution of Wealth'); 
label_x('time'); label_y('real wealth'); 
formataxis(gca); set(gca,'box','off');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end    
cc=cc+1;
delete([path_g '/*.eps']);
return
%% Diagnostics
% Plotting Residual Functions - Internal Use
figure;
plot(Periods, Z_M,'b'); title('Excess Money Supply'); drawnow;

figure;
plot(Periods, Z_Y,'b'); title('Excess Supply of Goods'); drawnow;

figure;
plot(Periods, Z_S,'b'); title('Excess Supply of Assets'); drawnow;

figure;
plot(Periods, p_t(Periods),'LineWidth',2,'Color',color1_ss); 
title('Price Path'); drawnow; 
orient landscape;

figure;
plot(Periods, Y_t,'b'); hold on;
plot(Periods, C_t,'r--');
legend('Output','Consumption'); 
