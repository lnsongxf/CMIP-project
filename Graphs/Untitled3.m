imprime      = @(x) print( gcf, '-depsc2', [path_g filesep x]);
imprpdf      = @(x) eps2pdf( [path_g filesep x '.eps']);
formataxis   = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 18, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
formatlegend = @(x) set(legend, 'Location', x, 'Orientation', 'Vertical', 'Box', 'On', 'Fontsize', 20, 'Fontangle', 'normal');
label_x      = @(x) xlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 14.25,'interpreter','latex');
label_y      = @(x) ylabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 14.25,'interpreter','latex');
label_z      = @(x) zlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 14.25,'interpreter','latex');
otitle       = @(x) title(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15,'interpreter','latex');


figure(cc);
plot(Periods(1:t_plotmax), rs_t(1:t_plotmax)*100,'LineWidth',4,'Color',color1_ss); hold on;
plot(Periods(1:t_plotmax), rb_t(1:t_plotmax)*100,'--o','LineWidth',4,'Color',color1_shock);  grid on;
plot(Periods(1:t_plotmax), rs_ss*100+0*rb_t(1:t_plotmax)*100,'k:','LineWidth',3); 
axis tight; [yout]=get(gca,'ylim');
line([T_pre+1 T_pre+1],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
line([T_post T_post],[yout(1) yout(2)],'Color','k','LineWidth',2,'LineStyle','--');
h=patch([T_pre+1 T_post T_post T_pre+1],[yout(1) yout(1) yout(2) yout(2)],[0.9 0.9 0.9]); alpha(h,0.2);
legend('r^a','r^l','r^{a}_{ss}');
label_x('time'); label_y('$\%$');
formataxis(gca); 
formatlegend('southeast');
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end