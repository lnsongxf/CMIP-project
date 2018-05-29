clear all; clc;

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


t_plotmax = 80;
Periods = (1:1:100);
color1_ss   =[0.1 0.1 0.6]; %color2_ss = color1_ss; color3_ss = color1_ss;
color2_ss   = [220 50 0]/255;%[0.8 0.4 0.4];
color3_ss   = [255 200 0]/255;%[0.6 0.8 0.6];

load('DDfp2.mat');
endopath1 = endopath;
steady1 = steady;
load('DDfp_ext.mat');
endopath2 = endopath;
steady2 = steady;
load('DDbp_ext.mat');
endopath3 = endopath;
steady3 = steady;

T_pre = 15       ; % Time Lapse prior to transition
T_post= 25       ; % Time Lapse after transition

figure(1);
hh=plot(Periods(1:t_plotmax)-1, [(endopath1.Y_t(1:t_plotmax)'/steady1.Y_ss-1) ...
                               (endopath2.Y_t(1:t_plotmax)'/steady2.Y_ss-1) ...
                                (endopath3.Y_t(1:t_plotmax)'/steady3.Y_ss-1)]*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',3,'Color',color1_ss);
set(hh(2),'LineWidth',3,'Color',color2_ss,'marker','o','linestyle','--');
set(hh(3),'LineWidth',3,'Color',color3_ss,'marker','x','linestyle',':','markersize',8)
axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
legend('wo/ADE','w/ADE','w/ADE + \Delta r');
label_x('time'); label_y('$\%$');
imprime(['comp_yt']);
imprpdf(['comp_yt']);

figure(2);
hh=plot(Periods(1:t_plotmax)-1, [(endopath1.B_t(1:t_plotmax)'/steady1.B_ss-1) ...
                                (endopath2.B_t(1:t_plotmax)'/steady2.B_ss-1) ...
                                (endopath3.B_t(1:t_plotmax)'/steady3.B_ss-1)]*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',3,'Color',color1_ss);
set(hh(2),'LineWidth',3,'Color',color2_ss,'marker','o','linestyle','--');
set(hh(3),'LineWidth',3,'Color',color3_ss,'marker','x','linestyle',':','markersize',8)
axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
% legend('wo/ADE','w/ADE','w/ADE + \Delta r');
label_x('time'); label_y('$\%$');
imprime(['comp_Bt']);
imprpdf(['comp_Bt']);

figure(3);
hh=plot(Periods(1:t_plotmax)-1, [endopath1.rs_t(1:t_plotmax)' ...
                               endopath2.rs_t(1:t_plotmax)' ...
                               endopath3.rs_t(1:t_plotmax)' ...
                               endopath3.rb_t(1:t_plotmax)']*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',3,'Color',color1_ss);
set(hh(2),'LineWidth',3,'Color',color2_ss,'marker','o','linestyle','--');
set(hh(3),'LineWidth',3,'Color',color3_ss,'marker','x','linestyle',':','markersize',8)
set(hh(4),'LineWidth',3,'Color',color3_ss,'marker','x','linestyle',':','markersize',8)

axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
% legend('wo/ADE','w/ADE','w/ADE + \Delta r');
label_x('time'); label_y('$\%$');
imprime(['comp_rt']);
imprpdf(['comp_rt']);


