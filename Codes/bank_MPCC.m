%% it is a shorcut for running banks commands in MPPC_generalpolicy.m
printit = 0;

i_m_ss    = -0.000;   %-0.005  ;
isp_ss    = 0.0512;   %0.02  ;

% Interbank Market block
varrho= 0.10 ; % Reserve Requirement
barlam= 2.5  ; % Efficiency interbank market
eta   = 0.5  ; % Bargaining
omega = 0.08 ; % Average size of shock - interbank

% Steady State Targets:
% For ploting Tests
imonmu_ss = varrho;
% i_m_ss    = 0.00;   %-0.005  ;
% isp_ss    = 0.03;   %0.02  ;


% For Plots
step       = 0.001;
surplus_cut=(varrho-omega*(1-varrho));
deficit_cut=(varrho+omega*(1-varrho));
delta_cut=deficit_cut-surplus_cut;
mu_vec=[(surplus_cut-0.1*delta_cut:step:surplus_cut)...
    (surplus_cut+step:step:deficit_cut-step)...
    (deficit_cut:step:deficit_cut+0.1*delta_cut)];

M_aux     = 100        ;
fedeqaux  = 000        ;

poldaux = isp_ss; % typical policy spread
r_m_aux = i_m_ss;


MPCC_bankblock_ii;

MPPC_report_shocks;

% Interbank Market Solutions
MPCC_interbank_vecs;

% Interbank market plots
plotib=0;
if plotib==1
    MPCC_interbank_plots;
end

MPCC_nominalimplementation;