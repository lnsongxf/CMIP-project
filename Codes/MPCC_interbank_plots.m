% Cutoffs - In liquidity space
printit=1;

% In money-multiplier space
mm_surplus_cut=(varrho-omega*(1-varrho))^(-1);
mm_deficit_cut=(varrho+omega*(1-varrho))^(-1);

%% Begin main figures - Liquidity Space
% Liquidity Ratio to Interbank Tightness
figure(cc);
plot(mu_vec,log(theta(mu_vec)),'LineWidth',4);
xlabel('$\Lambda$','interpreter','latex'); ylabel('log($\theta$)','interpreter','latex');
ftitle=title('Liquidiy Ratio and Interbank Tightness','interpreter','latex'); grid on;
axis tight;
[yout]=get(gca,'ylim');
line([surplus_cut surplus_cut],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([deficit_cut deficit_cut],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([varrho varrho],[yout(1) yout(2)],'Color',[0 0 0.2],'LineWidth',2);
% Labels
MPCC_autoplot;
text(varrho , yout(1),'$\varrho$','interpreter','latex','HorizontalAlignment','center','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
text(deficit_cut , yout(1),'$\varrho+\delta (1-\varrho)$','interpreter','latex','HorizontalAlignment','right','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
text(surplus_cut , yout(1),'$\varrho-\delta (1-\varrho)$','interpreter','latex','HorizontalAlignment','left','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
if printit==1
    orient landscape;
    saveas(gcf,'F_thetaLambda','pdf');
end
cc=cc+1;

% Final vs. Initial Tightness
figure(cc);
hold on;
grid on;
plot(log(theta_vec),log(bartheta_vec),'-','LineWidth',4);
plot(log(theta_vec),log(theta_vec),'k:','LineWidth',4);
ftitle=title('Final and Initial Interbank Tightness','interpreter','latex'); grid on;
axis tight;
% MPCC_autoplot;
xlabel('$log(\theta)$','interpreter','latex'); axis tight;
set( gca                       , ...
     'FontSize'   , 20 , ...
     'FontName'   , 'Helvetica' );
set(ftitle, ...
     'FontName'   , 'AvantGarde');
set( ftitle                    , ...
     'FontSize'   , 20          , ...
     'FontWeight' , 'bold',...      
     'Fontname'   , 'Courier');
h=legend('$log(\bar{\theta})$','$45^{o}$');
set(h,'Interpreter','latex','FontSize',20);
if printit==1   
    orient landscape;
    saveas(gcf,'F_thetatheta','pdf');
end
cc=cc+1;

figure(cc);
plot(log(theta_vec),chi_p_vec,'--','LineWidth',4); hold on;
plot(log(theta_vec),chi_m_vec,'LineWidth',4); grid on;
plot(log(theta_vec),isp_ss+0*chi_m_vec,':','LineWidth',4); hold on;
MPCC_autoplot;
ftitle=title('Tightness and Average Interbank Rates','interpreter','latex'); grid on;
axis tight;
% MPCC_autoplot;
xlabel('$log(\theta)$','interpreter','latex'); axis tight;
set( gca                       , ...
     'FontSize'   , 20 , ...
     'FontName'   , 'Helvetica' );
set(ftitle, ...
     'FontName'   , 'AvantGarde');
set( ftitle                    , ...
     'FontSize'   , 20          , ...
     'FontWeight' , 'bold',...      
     'Fontname'   , 'Courier');
h=legend('$\chi^{+}$','$\chi^{-}$','$\iota$');
set(h,'Interpreter','latex','FontSize',20);
if printit==1
    orient landscape;
    saveas(gcf,'F_thetachi','pdf');
end
cc=cc+1;

figure(cc);
hold on;
plot(log(theta_vec),r_d_vec,'LineWidth',4); grid on;
plot(log(theta_vec),r_b_vec,'--','LineWidth',4); hold on;
plot(log(theta_vec),r_b_vec*0+poldaux+r_m_aux,':','LineWidth',4); hold on;
plot(log(theta_vec),r_b_vec*0+r_m_aux,'-.','LineWidth',4); hold on;
xlabel('$log(\theta)$','interpreter','latex'); xlim([log(theta_vec(end)) log(theta_vec(1))]);
title('Equilibrium Rates');
MPCC_autoplot;
h=legend('$i^{dw}$','$i^{m}$','$i^{l}$','$i^{a}$');
set(h,'Interpreter','latex','FontSize',20);
if printit==1
    orient landscape;
    saveas(gcf,'F_ratestheta','pdf');
end
cc=cc+1;

%% Figures in terms of Liquidity Ratio
figure(cc);
hold on;
plot(mu_vec,r_d_vec,'LineWidth',4); grid on;
plot(mu_vec,r_b_vec,'--','LineWidth',4); 
plot(mu_vec,r_b_vec*0+poldaux+r_m_aux,':','LineWidth',4); 
plot(mu_vec,r_b_vec*0+r_m_aux,'-.','LineWidth',4); 
xlabel('$\Lambda$','interpreter','latex'); 
ftitle=title('Equilibrium Rates and Liquidity Ratio','interpreter','latex');
MPCC_autoplot;
h=legend('$i^{a}$','$i^{l}$','$i^{dw}$','$i^{m}$');
set(h,'Interpreter','latex','FontSize',20);
axis tight;
[yout]=get(gca,'ylim');
% ylim([yout(1)*1.1 yout(2)*1.05]);
line([surplus_cut surplus_cut],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([deficit_cut deficit_cut],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([varrho varrho],[yout(1) yout(2)],'Color',[0 0 0.2],'LineWidth',2);
text(varrho , yout(1),'$\varrho$','interpreter','latex','HorizontalAlignment','center','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
text(deficit_cut , yout(1),'$\varrho+\delta (1-\varrho)$','interpreter','latex','HorizontalAlignment','right','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
text(surplus_cut , yout(1),'$\varrho-\delta (1-\varrho)$','interpreter','latex','HorizontalAlignment','left','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
if printit==1
    orient landscape;
    saveas(gcf,'F_rateslambda','pdf');
end
cc=cc+1;

figure(cc);
plot(mu_vec,fedbal_p_vec,'LineWidth',4); hold on;
plot(mu_vec,fedbal_m_vec,'LineWidth',4); hold off; grid on;
ftitle=title('Surplus and Deficit Positions at CB','interpreter','latex')
xlabel('$\Lambda$','interpreter','latex');
MPCC_autoplot;
axis tight;
h=legend('$B^{+}$','$B^{-}$'); 
set(h,'Interpreter','latex','FontSize',20);
[yout]=get(gca,'ylim');
% ylim([yout(1)*0.9 yout(2)*1.05]);
line([surplus_cut surplus_cut],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([deficit_cut deficit_cut],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([varrho varrho],[yout(1) yout(2)],'Color',[0 0 0.2],'LineWidth',2);
text(varrho , yout(1),'$\varrho$','interpreter','latex','HorizontalAlignment','center','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
text(deficit_cut , yout(1),'$\varrho+\delta (1-\varrho)$','interpreter','latex','HorizontalAlignment','right','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
text(surplus_cut , yout(1),'$\varrho-\delta (1-\varrho)$','interpreter','latex','HorizontalAlignment','left','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
if printit==1
    orient landscape;
    saveas(gcf,'F_balancesLambda','pdf');
end
cc=cc+1;

figure(cc);
plot(mu_vec,dwprofs_vec,'LineWidth',4); hold on; grid on;
plot(mu_vec,bsprofs_vec,'LineWidth',4); 
plot(mu_vec,totprofs_vec,'LineWidth',4);
ftitle=title('Composition of CB profits over Real Savings','interpreter','latex')
xlabel('$\Lambda$','interpreter','latex');
MPCC_autoplot; axis tight;
h=legend('from discount loans','from OMO arbitrage','$\pi^{cb}$ (total)');
set(h,'Interpreter','latex','FontSize',20);
[yout]=get(gca,'ylim');
% ylim([yout(1)*0.9 yout(2)*1.05]);
line([surplus_cut surplus_cut],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([deficit_cut deficit_cut],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([varrho varrho],[yout(1) yout(2)],'Color',[0 0 0.2],'LineWidth',2);
text(varrho , yout(1),'$\varrho$','interpreter','latex','HorizontalAlignment','center','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
text(deficit_cut , yout(1),'$\varrho+\delta (1-\varrho)$','interpreter','latex','HorizontalAlignment','right','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
text(surplus_cut , yout(1),'$\varrho-\delta (1-\varrho)$','interpreter','latex','HorizontalAlignment','left','VerticalAlignment','cap','FontSize',20,'FontWeight','bold')
if printit==1
    orient landscape;
    saveas(gcf,'F_profitmarginLambda','pdf');
end
cc=cc+1;

figure(cc);
plot(mu_vec,dwprofs_vec.*mu_vec,'LineWidth',3); hold on; grid on;
plot(mu_vec,bsprofs_vec.*mu_vec,'LineWidth',3,'LineStyle',':'); hold on; grid on;
plot(mu_vec,totprofs_vec.*mu_vec,'LineWidth',3,'LineStyle','--'); hold on; grid on;
legend('Discount Window Profits','Balance Sheet Profits','Total');
[yout]=get(gca,'ylim');
line([surplus_cut surplus_cut],[yout(1) yout(2)]);
line([deficit_cut deficit_cut],[yout(1) yout(2)]);
line([varrho varrho],[yout(1) yout(2)]);
title('FED profits per unit of Money');
xlabel('Rates as Function of Liquidity Ratio \Lambda');
if printit==1
    orient landscape;
    saveas(gcf,'F_pi_moneyLambda','pdf');
end
cc=cc+1;

figure(cc);
plot(mu_vec,totprofs_vec,'LineWidth',3); hold on; grid on;
[yout]=get(gca,'ylim');
line([surplus_cut surplus_cut],[yout(1) yout(2)]);
line([deficit_cut deficit_cut],[yout(1) yout(2)]);
line([varrho varrho],[yout(1) yout(2)]);
title('Check: profits Total');
xlabel('Rates as Function of Liquidity Ratio \Lambda');
if printit==1
    orient landscape;
    saveas(gcf,'F_piLambda','pdf');
end

% Checks: Interbank Block
r_d_vecaux=r_b_vec.*(1-mu_vec)+r_m_aux.*mu_vec+dwprofs_vec;
figure(cc);
plot(mu_vec,r_b_vec.*(1-mu_vec),'-','LineWidth',3); hold on; 
plot(mu_vec,r_m_aux.*mu_vec,'-.','LineWidth',3); hold on;
plot(mu_vec,r_d_vec,':','LineWidth',3);
plot(mu_vec,dwprofs_vec,':','LineWidth',3);
[yout]=get(gca,'ylim');
line([surplus_cut surplus_cut],[yout(1) yout(2)]);
line([deficit_cut deficit_cut],[yout(1) yout(2)]);
line([varrho varrho],[yout(1) yout(2)]);
grid on; axis tight;
title('Check 1');

% plot(mu_vec,dwprofs_vec,':','LineWidth',3); grid on;

zeroprof_check=r_b_vec.*(1-mu_vec)+r_m_aux.*mu_vec-r_d_vec-dwprofs_vec;
plot(mu_vec,zeroprof_check,':','LineWidth',3); grid on;
cc=cc+1;

figure(cc);
fin_profits=r_b_vec.*(1-mu_vec)+r_m_aux.*mu_vec-r_d_vec; hold on;
plot(mu_vec,fin_profits,'-.','LineWidth',3); hold on;
plot(mu_vec,dwprofs_vec,'-','LineWidth',3); hold on; 
title('Discount-Window Profits vs. Bank Portfolio Return');
legend('Discount-Window Profits','Bank Portfolio Return') ;
[yout]=get(gca,'ylim');
line([surplus_cut surplus_cut],[yout(1) yout(2)]);
line([deficit_cut deficit_cut],[yout(1) yout(2)]);
line([varrho varrho],[yout(1) yout(2)]);
grid on; axis tight;
title('Check 2');
% plot(mu_vec,r_d_vec,':','LineWidth',3);
cc=cc+1;

%% Dealing with the Zero-Lower Bound on Deposit Rates
%% Figures in terms of Liquidity Ratio
zlb_index    =find((r_d_vec<=0),1,'first');
r_d_vec_aux = r_d_vec;
r_b_vec_aux = r_b_vec;
dwprofs_vec_aux=dwprofs_vec  ;
bsprofs_vec_aux=bsprofs_vec  ;
totprofs_vec_aux=totprofs_vec;
if zlb_index<length(mu_vec);
    r_d_vec_aux(zlb_index:end)=r_d_vec_aux(zlb_index);
    r_b_vec_aux(zlb_index:end)=r_b_vec_aux(zlb_index);
    dwprofs_vec_aux(zlb_index:end)=dwprofs_vec(zlb_index)  ;
    bsprofs_vec_aux(zlb_index:end)=r_b_vec_aux(zlb_index)*mu_vec(zlb_index:end)-r_m_aux*mu_vec(zlb_index);
    totprofs_vec_aux(zlb_index:end)=bsprofs_vec_aux(zlb_index:end)+dwprofs_vec_aux(zlb_index:end);
end

figure(cc);
plot(mu_vec,r_b_vec_aux*0+poldaux+r_m_aux,'k.-'); hold on;
plot(mu_vec,r_b_vec*0+r_m_aux,'k+'); hold on;
plot(mu_vec,r_b_vec_aux,':','LineWidth',3); hold on;
plot(mu_vec,r_d_vec_aux,'LineWidth',3); grid on; legend('DW','RoR','Borrowing','Lending Rate');
xlabel('Liquidity Ratio \Lambda'); 
title('Borrowing\Lending Rates');
axis tight;
[yout]=get(gca,'ylim');
ylim([yout(1)*0.9 yout(2)*1.05]);
line([surplus_cut surplus_cut],[yout(1) yout(2)],'Color','k','LineWidth',2);
line([deficit_cut deficit_cut],[yout(1) yout(2)],'Color','k','LineWidth',2);
line([mu_vec(zlb_index) mu_vec(zlb_index)],[yout(1) yout(2)],'Color','k','LineWidth',2);
line([varrho varrho],[yout(1) yout(2)],'Color',[0 0 0.2],'LineWidth',2);
if printit==1
    orient landscape;
    saveas(gcf,'F_rateslambda_zlb','pdf');
end
cc=cc+1;

figure(cc);
plot(mu_vec,dwprofs_vec_aux,'LineWidth',3); hold on; grid on;
plot(mu_vec,bsprofs_vec_aux,'LineWidth',3); 
plot(mu_vec,totprofs_vec_aux,'LineWidth',3); 
legend('Discount Window Profits','Balance Sheet Profits','Total CB Operation Profits');
[yout]=get(gca,'ylim');
line([surplus_cut surplus_cut],[yout(1) yout(2)],'Color','k','LineWidth',2);
line([deficit_cut deficit_cut],[yout(1) yout(2)],'Color','k','LineWidth',2);
line([varrho varrho],[yout(1) yout(2)],'Color',[0 0 0.2],'LineWidth',2);
% line([mu_vec(zlb_index) mu_vec(zlb_index)],[yout(1) yout(2)],'Color','k','LineWidth',2);
title('CB profits/Private Savings');
xlabel('\Lambda'); axis tight
pseries={'dwprofs_vec_aux'}     ;
MPCC_autoplot                   ;
cc=cc+1;
