%% Computing Inflation and Monetary Policy Variables
% Compute Indices where it doesn't matter
inf_target = 0.00            ;   % Inflation Target
i_target   = rs_ss+inf_target;   % Nominal Rate Target
 
% Construct the index of admissible values

r_m_aux = i_m_ss;
MPCC_bankblock_ii;
MPCC_interbank_vecs;

mu_index=(mu_vec>surplus_cut & mu_vec<deficit_cut & mu_vec<mu_vec(index0));
[~,index_aux]=find(mu_index==1,1,'first');

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Policy Spread
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MB_ss1= mu_ss_target*D_ss          ; 
MB_t  = MB_ss1*ones(T+1,1)         ;  % constant along transition
ND_t  = D_t                        ;
mub_t = mu_ss_target*ones(T+1,1)   ;
mu_t  = mu_ss_target*ones(T+1,1)   ;
P_t   = ones(T+1,1)                ;
dP_t  = inf_target*ones(T+1,1)     ;
id_t  = i_target*ones(T+1,1)       ; 

for tt=1:T-1    
    if list(1)==1 
        r_m_aux = rr_t(tt+1);
        MPCC_bankblock_ii;
        MPCC_interbank_vecs;
        Dr=r_b_vec - r_d_vec;
        index0=zlb_index(tt+1);
        mu_index=(mu_vec>surplus_cut&mu_vec<deficit_cut & mu_vec<mu_vec(index0));
    end
    if rsp_t(tt+1) > Dr(index0)
        mub_t(tt+1) = interp1(Dr(mu_index),mu_vec(mu_index),rsp_t(tt+1),'pchip');
        id_t(tt+1)  = interp1(mu_vec(mu_index), r_d_vec(mu_index) ,mub_t(tt+1),'pchip');
    else
        mub_t(tt+1) = mu_vec(index0);
        id_t(tt+1)  = 0;
    end
    dP_t(tt+1)  = id_t(tt+1)-rs_t(tt);
    P_t(tt+1)   = P_t(tt)*(1+dP_t(tt+1)*dt);
    ND_t(tt+1)  = D_t(tt)*P_t(tt+1);
    MB_t(tt+1)  = mub_t(tt+1)*D_t(tt);
end
% M_t   = [D_ss;D_t'].*(mub_t +0.1*ones(T+1,1).*(mub_t==[mu_vec(zlb_index)';mu_vec(zlb_index(end))]));
% M0_t  = M_t-MB_t;

if list(1)==1 % constant level of money
    M_t  = MB_ss1*ones(T+1,1)         ;  % constant along transition
    M0_t = M_t - MB_t                 ;
else
    M_t   = [D_ss;D_t'].*(mub_t +0.1*ones(T+1,1).*(mub_t==mu_vec(index0)));
    M0_t  = M_t-MB_t;
end

index_1=(2:T+1);

figure(cc);
plot(index_1-2,M0_t(index_1-1)-M0_t(1),'LineWidth',4,'color',color1_ss); hold on;
plot(index_1-2,MB_t(index_1-1)-MB_t(1),'--o','LineWidth',4,'color',color1_shock); grid on;
xlim([0 t_plotmax]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
% otitle('Cash Holdings');
label_x('time');
legend('Currency','Reserves'); formatlegend('Best');axis tight;
formataxis(gca);  
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;


figure(cc);
plot(index_1-2,(P_t(index_1-1)),'LineWidth',4,'color',color1_ss);  grid on;
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
plot(index_1-2,id_t(index_1-1)*100,'color',[220 85 0]/255,'LineStyle','-','LineWidth',4); 
plot(index_1-2,rs_t(index_1-1)*100,'color',[255 195 0]/255,'LineWidth',4);
plot(index_1-2,inf_target*100+0*rs_t(index_1-1),'k-.','LineWidth',2); grid on;
ylim([-1 2])
[yout]=get(gca,'ylim');
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',1,'LineStyle','-');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
% otitle('Fisher Equation Decomp.');
legend('Inflation','Nominal Rate','Real Rate','Inflation Target'); 
label_x('time'); label_y('$\%$');
formataxis(gca); xlim([0 t_plotmax]);

formatlegend('northeast');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;
