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
color1_ss   =[0.1 0.1 0.8]; %color2_ss = color1_ss; color3_ss = color1_ss;
color2_ss   = [220 50 0]/255;%[0.8 0.4 0.4];
color3_ss   = [200 200 200]/255;%[0.6 0.8 0.6];
color4_ss   = [200 255 0]/255;%[0.6 0.8 0.6];
color5_ss   = [255 255 255]*0;%/255;%[0.6 0.8 0.6];

load('DDcc.mat');
s_plotmax = round(numel(s_vec)/4);
[X,Y]     = meshgrid((1:t_plotmax)',s_vec(1:s_plotmax));
ds        = s_vec(2)-s_vec(1);

endopath1 = endopath;
steady1   = steady;
Es1       = sum(endopath.f_t(:,1:t_plotmax).*repmat(s_vec,1,t_plotmax))*ds;
Vs1       = sum(endopath.f_t(:,1:t_plotmax).*(bsxfun(@minus,s_vec,Es1).^2))*ds;

load('DDccrsp.mat');
endopath2 = endopath;
steady2 = steady;
Es2       = sum(endopath.f_t(:,1:t_plotmax).*repmat(s_vec,1,t_plotmax))*ds;
Vs2       = sum(endopath.f_t(:,1:t_plotmax).*(bsxfun(@minus,s_vec,Es2).^2))*ds;

load('DDccrspendo.mat');
endopath3 = endopath;
steady3   = steady;
Es3       = sum(endopath.f_t(:,1:t_plotmax).*repmat(s_vec,1,t_plotmax))*ds;
Vs3       = sum(endopath.f_t(:,1:t_plotmax).*(bsxfun(@minus,s_vec,Es3).^2))*ds;

load('DDcc_rss.mat');
endopath4 = endopath;
steady4   = steady;
Es4       = sum(endopath.f_t(:,1:t_plotmax).*repmat(s_vec,1,t_plotmax))*ds;
Vs4       = sum(endopath.f_t(:,1:t_plotmax).*(bsxfun(@minus,s_vec,Es4).^2))*ds;

load('DDccrsp_rss.mat');
endopath5 = endopath;
steady5   = steady;
Es5       = sum(endopath.f_t(:,1:t_plotmax).*repmat(s_vec,1,t_plotmax))*ds;
Vs5       = sum(endopath.f_t(:,1:t_plotmax).*(bsxfun(@minus,s_vec,Es5).^2))*ds;

load('DDccrspendo_rss.mat');
endopath6 = endopath;
steady6   = steady;
Es6       = sum(endopath.f_t(:,1:t_plotmax).*repmat(s_vec,1,t_plotmax))*ds;
Vs6       = sum(endopath.f_t(:,1:t_plotmax).*(bsxfun(@minus,s_vec,Es6).^2))*ds;

T_pre = 1       ; % Time Lapse prior to transition
T_post= 16       ; % Time Lapse after transition

figure(1);
hh=plot(Periods(1:t_plotmax)-1, [(endopath1.Y_t(1:t_plotmax)'/steady1.Y_ss-1) ...
                               (endopath2.Y_t(1:t_plotmax)'/steady2.Y_ss-1) ...
                                (endopath3.Y_t(1:t_plotmax)'/steady3.Y_ss-1) ...
                                (endopath4.Y_t(1:t_plotmax)'/steady4.Y_ss-1) ...
                                (endopath5.Y_t(1:t_plotmax)'/steady5.Y_ss-1) ...
                                (endopath6.Y_t(1:t_plotmax)'/steady6.Y_ss-1) ]*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',1.5,'Color',color1_ss);
set(hh(2),'LineWidth',1.5,'Color',color2_ss);
set(hh(3),'LineWidth',1.5,'Color',color5_ss);
set(hh(4),'LineWidth',1,'Color',color1_ss,'marker','o','linestyle','--')
set(hh(5),'LineWidth',1,'Color',color2_ss,'marker','o','linestyle','--')
set(hh(6),'LineWidth',1,'Color',color5_ss,'marker','o','linestyle','--','markersize',6)
% set(hh(7),'LineWidth',1.5,'Color','g');%,'marker','s','linestyle',':','markersize',8)
% set(hh(8),'LineWidth',1,'Color','g','marker','o','linestyle','--','markersize',6)


axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
legend('CC','CC+\Delta r','CC+\Delta r endo','CC\_rss','CC\_rss+\Delta r','CC\_rss+\Delta r endo');
label_x('time'); label_y('$\%$');
imprime(['comp_cc_yt']);
imprpdf(['comp_cc_yt']);

%%
figure(2);
hh=plot(Periods(1:t_plotmax)-1, [(endopath1.B_t(1:t_plotmax)'/steady1.B_ss-1) ...
                                (endopath2.B_t(1:t_plotmax)'/steady2.B_ss-1) ...
                                (endopath3.B_t(1:t_plotmax)'/steady3.B_ss-1) ...
                                (endopath4.B_t(1:t_plotmax)'/steady4.B_ss-1) ...
                                (endopath5.B_t(1:t_plotmax)'/steady5.B_ss-1) ...
                                (endopath6.B_t(1:t_plotmax)'/steady6.B_ss-1) ]*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',1.5,'Color',color1_ss);
set(hh(2),'LineWidth',1.5,'Color',color2_ss);
set(hh(3),'LineWidth',1,'Color',color1_ss,'marker','o','linestyle','--')
set(hh(4),'LineWidth',1,'Color',color2_ss,'marker','o','linestyle','--')
set(hh(5),'LineWidth',1.5,'Color',color5_ss);%,'marker','x','linestyle',':','markersize',8)
set(hh(6),'LineWidth',1,'Color',color5_ss,'marker','o','linestyle','--','markersize',6)
% set(hh(7),'LineWidth',1.5,'Color','g');%,'marker','s','linestyle',':','markersize',8)
% set(hh(8),'LineWidth',1,'Color','g','marker','o','linestyle','--','markersize',6)

axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
% legend('wo/ADE','w/ADE','w/ADE + \Delta r');
label_x('time'); label_y('$\%$');
imprime(['comp_cc_Bt']);
imprpdf(['comp_cc_Bt']);


figure(3);
hh=plot(Periods(1:t_plotmax)-1, [(endopath1.D_t(1:t_plotmax)'/steady1.D_ss-1) ...
                                (endopath2.D_t(1:t_plotmax)'/steady2.D_ss-1) ...
                                (endopath3.D_t(1:t_plotmax)'/steady3.D_ss-1) ...
                                (endopath4.D_t(1:t_plotmax)'/steady4.D_ss-1) ...
                                 (endopath5.D_t(1:t_plotmax)'/steady5.D_ss-1) ...
                                 (endopath6.D_t(1:t_plotmax)'/steady6.D_ss-1)]*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',1.5,'Color',color1_ss);
set(hh(2),'LineWidth',1.5,'Color',color2_ss);
set(hh(3),'LineWidth',1.5,'Color',color5_ss);
set(hh(4),'LineWidth',1,'Color',color1_ss,'marker','o','linestyle','--')
set(hh(5),'LineWidth',1,'Color',color2_ss,'marker','o','linestyle','--')
set(hh(6),'LineWidth',1,'Color',color5_ss,'marker','o','linestyle','--','markersize',6)
% set(hh(7),'LineWidth',1.5,'Color','g');%,'marker','s','linestyle',':','markersize',8)
% set(hh(8),'LineWidth',1,'Color','g','marker','o','linestyle','--','markersize',6)

axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
% legend('wo/ADE','w/ADE','w/ADE + \Delta r');
label_x('time'); label_y('$\%$');
imprime(['comp_cc_Dt']);
imprpdf(['comp_cc_Dt']);

figure(4);
hh=plot(Periods(1:t_plotmax)-1, [endopath1.rb_t(1:t_plotmax)'-endopath1.rs_t(1:t_plotmax)' ...
                               endopath2.rb_t(1:t_plotmax)' - endopath2.rs_t(1:t_plotmax)' ...
                               endopath3.rb_t(1:t_plotmax)' - endopath3.rs_t(1:t_plotmax)' ...
                               endopath4.rb_t(1:t_plotmax)' - endopath4.rs_t(1:t_plotmax)' ...
                               endopath5.rb_t(1:t_plotmax)' - endopath5.rs_t(1:t_plotmax)' ...
                               endopath6.rb_t(1:t_plotmax)' - endopath6.rs_t(1:t_plotmax)']*100); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',1.5,'Color',color1_ss);
set(hh(2),'LineWidth',1.5,'Color',color2_ss);
set(hh(3),'LineWidth',1.5,'Color',color5_ss);
set(hh(4),'LineWidth',1,'Color',color1_ss,'marker','o','linestyle','--')
set(hh(5),'LineWidth',1,'Color',color2_ss,'marker','o','linestyle','--')
set(hh(6),'LineWidth',1,'Color',color5_ss,'marker','o','linestyle','--','markersize',6)
% set(hh(7),'LineWidth',1.5,'Color','g');%,'marker','s','linestyle',':','markersize',8)
% set(hh(8),'LineWidth',1,'Color','g','marker','o','linestyle','--','markersize',6)

axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
% legend('wo/ADE','w/ADE','w/ADE + \Delta r');
label_x('time'); label_y('$\%$');
imprime(['comp_cc_rspt']);
imprpdf(['comp_cc_rspt']);

%%
figure(5);
hh=plot(Periods(1:t_plotmax)-1, [Vs1' Vs2' Vs3' Vs4' Vs5' Vs6']); hold on;
plot(Periods(1:t_plotmax)-1, 0*endopath1.rb_t(1:t_plotmax),'k:','LineWidth',1); grid on;
set(hh(1),'LineWidth',1.5,'Color',color1_ss);
set(hh(2),'LineWidth',1.5,'Color',color2_ss);
set(hh(3),'LineWidth',1.5,'Color',color5_ss);
set(hh(4),'LineWidth',1,'Color',color1_ss,'marker','o','linestyle','--')
set(hh(5),'LineWidth',1,'Color',color2_ss,'marker','o','linestyle','--')
set(hh(6),'LineWidth',1,'Color',color5_ss,'marker','o','linestyle','--','markersize',6)
% set(hh(7),'LineWidth',1.5,'Color','g');%,'marker','s','linestyle',':','markersize',8)
% set(hh(8),'LineWidth',1,'Color','g','marker','o','linestyle','--','markersize',6)

axis tight; [yout]=get(gca,'ylim'); 
ylim([yout(1) yout(2)*1.1]);
[yout]=get(gca,'ylim'); 
line([T_pre T_pre],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post-1 T_post-1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre T_post-1 T_post-1 T_pre],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
formataxis(gca); %otitle('$Y_{t}$ deviation from $Y_{ss}$'); 
% legend('wo/ADE','w/ADE','w/ADE + \Delta r');
label_x('time'); label_y('$\%$');
imprime(['comp_cc_Vs']);
imprpdf(['comp_cc_Vs']);

delete(['\Graphs\*.eps'])
return
%%
[steady1.Y_ss steady2.Y_ss steady3.Y_ss steady4.Y_ss steady5.Y_ss steady6.Y_ss]
