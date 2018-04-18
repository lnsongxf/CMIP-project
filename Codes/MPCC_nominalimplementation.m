%% Computing Inflation and Monetary Policy Variables
% Compute Indices where it doesn't matter
inf_target = 0.00     ;          % Inflation Target
i_target   = rs_ss+inf_target;   % Nominal Rate Target
 
% Construct the index of admissible values
mu_index=(mu_vec>surplus_cut&mu_vec<deficit_cut);
[~,index_aux]=find(mu_index==1,1,'first');
if index_aux>1
    mu_index(index_aux-1)=1;
end
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

% Plot Map
figure
plot(mu_vec(mu_index),[r_b_vec(mu_index)' r_d_vec(mu_index)'],'LineWidth',3); hold on;
scatter(mu_ss_target,i_target,40,'MarkerEdgeColor','r','MarkerFaceColor','r');
ylabel('Nominal Rate i^{d}'); xlabel('Liquidity Ratio \Lambda'); grid on; hold on;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Monetary Rate Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MB_ss1= mu_ss_target*D_ss    ;
MB_t  = MB_ss1*ones(T+1,1)       ;  % constant along transition
ND_t  = D_t                      ;
mu_t  = mu_ss_target*ones(T,1)   ;
P_t   = ones(T+1,1)              ;
dP_t  = inf_target*ones(T,1)     ;
id_t  = i_target*ones(T,1)       ;

for tt=1:T-1
    dP_t(tt)=id_t(tt)-rb_t(tt);
    P_t(tt+1)=P_t(tt)*(1+dP_t(tt)*dt);
    ND_t(tt+1)=D_t(tt+1)*P_t(tt+1);
    mu_t(tt+1)=MB_t(tt+1)/ND_t(tt+1);
    id_t(tt+1)=max(interp1(mu_vec,r_d_vec,mu_t(tt+1),'linear'),0);
end
index_1=(1:T-1);

figure(cc);
plot(index_1,(P_t(index_1)),'LineWidth',3,'color',[.1 .1 .7]); title('Price Path - Constant Monetary Rule'); 
xlabel('Time','FontSize',25); ylabel('P_t','FontSize',25); title('Price Path'); grid on;
[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post T_post],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post T_post T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
ftitle=title('Price Path','interpreter','latex','fontsize',20); grid on;
axis tight;
formataxis(gca);
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc);
plot(index_1,dP_t(index_1)*100,'LineStyle',':','LineWidth',4,'color',[.1 .1 .7]); hold on;
plot(index_1,id_t(index_1)*100,'color',[220 85 0]/255,'LineStyle','--','LineWidth',4); 
plot(index_1,rb_t(index_1)*100,'color',[255 195 0]/255,'LineWidth',4); grid on;
plot(index_1,inf_target*100+0*rb_t(index_1),'k-.','LineWidth',2); grid on;
title('Components for Fisher Equation')
legend('Inflation','Nominal Rate','Real Rate','Inflation Target'); axis tight;
xlabel('Time','FontSize',25); ylabel('$\%$','interpreter','latex','FontSize',25);
[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post T_post],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post T_post T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
ftitle=title('Fisher Equation Decomposition','interpreter','latex','fontsize',20); grid on;
axis tight;
formataxis(gca);
formatlegend('Best');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inflation Target Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psipi=0.0355;
psiyi=0.0125;

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
    dP_t(tt)=id_t(tt)-rb_t(tt+1)                  ;
    P_t(tt+1)=P_t(tt)*(1+dP_t(tt)*dt)             ;
    ND_t(tt+1)=D_t(tt+1)*P_t(tt+1)                ;
    id_t(tt+1)=max(0,i_target+psipi*(dP_t(tt)-inf_target)+psiyi*(Y_t(tt)-Y_ss));    % Zero lower bound
    mu_t(tt+1)=interp1(r_d_aux,mu_aux,id_t(tt+1),'pchip');
end
index_1=(1:T-1);

figure(cc);
plot(index_1,(P_t(index_1)),'LineWidth',4,'color',[.1 .1 .7]); title('Price Path - Constant Monetary Rule'); 
xlabel('t','interpreter','latex','FontSize',25); ylabel('P_t','FontSize',25);%ylabel('$\%$','interpreter','latex','FontSize',25);
ftitle=title('Price Path','interpreter','latex','fontsize',20); grid on;
[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post T_post],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post T_post T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
axis tight;
formataxis(gca);

if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

figure(cc);
plot(index_1,inf_target*100+0*rb_t(index_1),'k-.','LineWidth',2);  hold on;
plot(index_1,dP_t(index_1+1)*100,'LineStyle',':','LineWidth',4);  hold on;
plot(index_1,id_t(index_1+1)*100,'color',[220 85 0]/255,'LineStyle','--','LineWidth',4); 
plot(index_1,rb_t(index_1+1)*100,'color',[255 195 0]/255,'LineWidth',4); grid on;
grid on;
xlabel('t','interpreter','latex','FontSize',25); ylabel('$\%$','interpreter','latex','FontSize',25);
ftitle=title('Fisher Equation Decomposition','interpreter','latex','fontsize',20); grid on;
axis tight;
[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post T_post],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post T_post T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]);
alpha(h,0.2);

h=legend('Inflation Target','Inflation','Nominal Rate','Real Rate'); 
formataxis(gca);axis tight;
formatlegend('Best');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;