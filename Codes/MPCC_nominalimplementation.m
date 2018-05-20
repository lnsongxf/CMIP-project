%% Computing Inflation and Monetary Policy Variables

% Interbank Market block 
varrho  = 0.5 ; % Reserve Requirement
barlam  = 2.5 ; % Efficiency interbank market
eta     = 0.5 ; % Bargaining
omega   = 0.42 ; % Average size of shock - interbank

% Steady State Targets:
imonmu_ss = varrho;
i_m_ss    = 0.00;
isp_ss    = 0.06;   

% For Plots
step       = 0.0001;
surplus_cut=(varrho-omega*(1-varrho));
deficit_cut=(varrho+omega*(1-varrho));
delta_cut=deficit_cut-surplus_cut;
mu_vec=[(surplus_cut-0.0001*delta_cut:step:surplus_cut)...
    (surplus_cut+step:step:deficit_cut-step)...
    (deficit_cut:step:deficit_cut+0.01*delta_cut)];

M_aux     = 100        ;
fedeqaux  = 000        ;

poldaux = isp_ss; % typical policy spread
r_m_aux = i_m_ss;   % this will be adjusted 


MPCC_bankblock_ii;
MPCC_interbank_vecs;
Dr  = r_b_vec-r_d_vec;

% Compute Indices where it doesn't matter
inf_target = 0.00     ;          % Inflation Target
i_target   = rs_ss+inf_target;   % Nominal Rate Target
 
% Construct the index of admissible values
mu_index=(mu_vec>surplus_cut&mu_vec<deficit_cut);
[~,index_aux]=find(mu_index==1,1,'first');
% if index_aux>1
%     mu_index(index_aux-1)=1;
% end
[~,index_aux]=find(mu_index==1,1,'last');
if index_aux<length(mu_index);
    mu_index(index_aux+1)=1;
end

% Find Policy
if i_target<=max(r_d_vec(mu_index)) & i_target>=min(r_d_vec(mu_index))
   mu_ss_target  = interp1(r_d_vec(mu_index),mu_vec(mu_index),i_target,'pchip');
   display('Liquidity ratio target is:'); 
   mu_ss_target
else
   display('Cannot implement that target'); 
end


figure(cc)
plot(mu_vec(mu_index),[r_b_vec(mu_index)' r_d_vec(mu_index)'],'LineWidth',3); hold on;
scatter(mu_ss_target,rs,40,'MarkerEdgeColor','r','MarkerFaceColor','r');
label_y('nominal rate'); label_x('Liquidity Ratio $\Lambda$'); grid on; hold on;
legend('i^{b}','i^{d}');
axis tight;
formataxis(gca); formatlegend('best');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Monetary Rate Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MB_ss1= mu_ss_target*D_ss    ;
MB_t  = MB_ss1*ones(T,1); %*cumprod(ones(T,1)-0.01)       ;  % constant along transition
ND_t  = D_t                      ;
mu_t  = mu_ss_target*ones(T,1)   ;
P_t   = ones(T+1,1)              ;
dP_t  = inf_target*ones(T,1)     ;
id_t  = i_target*ones(T,1)       ;

for tt=1:T-1
    dP_t(tt)=id_t(tt)-rs_t(tt);
    P_t(tt+1)=P_t(tt)*(1+dP_t(tt)*dt);
    ND_t(tt+1)=D_t(tt+1)*P_t(tt+1);
    mu_t(tt+1)=MB_t(tt+1)/ND_t(tt+1);
    id_t(tt+1)=max(interp1(mu_vec,r_d_vec,mu_t(tt+1),'linear'),0);
end
index_1=(2:T+1);

figure(cc);
plot(index_1-2,(P_t(index_1-1)),'LineWidth',3,'color',[.1 .1 .7]);  grid on;
[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
label_x('time'); label_y('Price index'); %otitle('Price Index Path');
formataxis(gca); axis tight; xlim([0 t_plotmax]);
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc);
plot(index_1-2,dP_t(index_1-1)*100,'LineStyle',':','LineWidth',4,'color',[.1 .1 .7]); hold on;
plot(index_1-2,id_t(index_1-1)*100,'color',[220 85 0]/255,'LineStyle','--','LineWidth',4); 
plot(index_1-2,rs_t(index_1-1)*100,'color',[255 195 0]/255,'LineWidth',4);
plot(index_1-2,inf_target*100+0*rs_t(index_1-1),'k-.','LineWidth',2); grid on;

[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
% otitle('Fisher Equation Decomp.');
legend('Inflation','Nominal Rate','Real Rate','Inflation Target'); axis tight;
label_x('time'); label_y('$\%$');
formataxis(gca); axis tight; xlim([0 t_plotmax]);

formatlegend('northeast');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inflation Target Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psipi=0.0125;
psiyi=0.4;

% Construc auxiliaries
index_rd=(r_b_vec>min(r_b_vec))& (r_b_vec<max(r_b_vec));
r_d_aux=r_d_vec(index_rd);
mu_aux=mu_vec(index_rd);

% Constructing time zero rate
MB_ss1= mu_ss_target*D_ss        ;
MB_t  = MB_ss1*ones(T+1,1)       ;
ND_t  = D_t                      ;
mu_t  = mu_ss_target*ones(T,1) ;
P_t   = ones(T+1,1)              ;
dP_t  = inf_target*ones(T,1)     ;
id_t  = i_target*ones(T,1)       ;

for tt=1:T-1
    dP_t(tt)=id_t(tt)-rs_t(tt)                  ;
    P_t(tt+1)=P_t(tt)*(1+dP_t(tt)*dt)             ;
    ND_t(tt+1)=D_t(tt+1)*P_t(tt+1)                ;
    id_t(tt+1)=max(0,i_target+psipi*(dP_t(tt)-inf_target)+psiyi*(Y_t(tt)-Y_ss)/Y_ss);    % Zero lower bound
    mu_t(tt+1)=interp1(r_d_aux,mu_aux,id_t(tt+1),'pchip');
end
index_1=(2:T+1);

figure(cc);
plot(index_1-2,(P_t(index_1-1)),'LineWidth',3,'color',[.1 .1 .7]);  grid on;
[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
label_x('time'); label_y('Price index'); %otitle('Price Index Path');
formataxis(gca); axis tight; xlim([0 t_plotmax]);
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc);
plot(index_1-2,dP_t(index_1-1)*100,'LineStyle',':','LineWidth',4,'color',[.1 .1 .7]); hold on;
plot(index_1-2,id_t(index_1-1)*100,'color',[220 85 0]/255,'LineStyle','--','LineWidth',4); 
plot(index_1-2,rs_t(index_1-1)*100,'color',[255 195 0]/255,'LineWidth',4);
plot(index_1-2,inf_target*100+0*rs_t(index_1-1),'k-.','LineWidth',2); grid on;

[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
% otitle('Fisher Equation Decomp.');
legend('Inflation','Nominal Rate','Real Rate','Inflation Target'); axis tight;
label_x('time'); label_y('$\%$');
formataxis(gca); axis tight; xlim([0 t_plotmax]);

formatlegend('northeast');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;