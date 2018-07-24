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


t_plotmax   = 80;
Periods     = (1:1:100);
color1_ss   =[0.1 0.1 0.8]; %color2_ss = color1_ss; color3_ss = color1_ss;
color2_ss   = [220 50 0]/255;%[0.8 0.4 0.4];
color3_ss   = [255 200 0]/255;%[0.6 0.8 0.6];
color4_ss   = [200 255 0]/255;%[0.6 0.8 0.6];

load('DDbp_zlb0.mat');
ds = s_vec(2)-s_vec(1);
fc_t = zeros(t_plotmax+1,1);  
for t=1:t_plotmax+1
    temp    = paths.s_bl_index_t(t); 
    fc_t(t) = sum(f_t(1:temp,t)*ds);
end
endopath.fc_t = fc_t;
endopath1     = endopath;
steady1       = steady;

load('DDbp_zlb2.mat');
ds = s_vec(2)-s_vec(1);
fc_t = zeros(t_plotmax+1,1);  
for t=1:t_plotmax+1
    temp    = paths.s_bl_index_t(t); 
    fc_t(t) = sum(f_t(1:temp,t)*ds);
%     C_t(t)  = f_t(:,t).*
end
endopath.fc_t = fc_t;
endopath2 = endopath;
steady2 = steady;

T_pre = 15        ; % Time Lapse prior to transition
T_post= 30       ; % Time Lapse after transition

%%
figure(1);
hh=plot(Periods(1:t_plotmax)-1, [(endopath1.Y_t(1:t_plotmax)'/steady1.C_ss-1) ...
                               (endopath2.Y_t(1:t_plotmax)'/steady2.C_ss-1)]*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',3,'Color',color1_ss);
set(hh(2),'LineWidth',3,'Color',color2_ss,'marker','o','linestyle','--');
axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
legend('no MP intervention','MP intervention'); formatlegend('best');
label_x('time'); label_y('$\%$');
imprime(['com_bp_yt']);
imprpdf(['comp_bp_yt']);


figure(2);
hh=plot(Periods(1:t_plotmax)-1, [endopath1.fc_t(1:t_plotmax) ...
                                 endopath2.fc_t(1:t_plotmax)]*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',3,'Color',color1_ss);
set(hh(2),'LineWidth',3,'Color',color2_ss,'marker','o','linestyle','--');
axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
legend('no MP intervention','MP intervention'); formatlegend('best');
label_x('time'); label_y('$\%$');
imprime(['com_bp_fct']);
imprpdf(['comp_bp_fct']);


%%
figure(3);
hh=plot(Periods(1:t_plotmax)-1, [(endopath1.B_t(1:t_plotmax)'/steady1.B_ss-1) ...
                                (endopath2.B_t(1:t_plotmax)'/steady2.B_ss-1)]*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',3,'Color',color1_ss);
set(hh(2),'LineWidth',3,'Color',color2_ss,'marker','o','linestyle','--');
axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca);
legend('no MP intervention','MP intervention'); formatlegend('best');
label_x('time'); label_y('$\%$');
imprime(['comp_bp_Bt']);
imprpdf(['comp_bp_Bt']);

figure(4);
hh=plot(Periods(1:t_plotmax)-1, [-endopath1.rs_t(1:t_plotmax)'+endopath1.rb_t(1:t_plotmax)' ...
                                -endopath2.rs_t(1:t_plotmax)'+endopath2.rb_t(1:t_plotmax)']*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',3,'Color',color1_ss);
set(hh(2),'LineWidth',3,'Color',color2_ss,'marker','o','linestyle','--');
axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca);
legend('no MP intervention','MP intervention'); formatlegend('best');
label_x('time'); label_y('$\%$');
imprime(['comp_bp_rt']);
imprpdf(['comp_bp_rt']);


