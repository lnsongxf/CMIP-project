%% Finding monetary policy instruments
% Interbank Market Equilibrium
MPCC_bankblock_ii;
MPCC_interbank_vecs;
Dr  = r_b_vec-r_d_vec;

mu_index=(mu_vec>surplus_cut & mu_vec<deficit_cut);
[~,index_aux] =find(mu_index==1,1,'first');

if index_aux<length(mu_index);
    mu_index(index_aux+1)=1;
end

mu_ss_target    = interp1(Dr(mu_index),mu_vec(mu_index),rsp_ss,'pchip');
[~,index_aux2]  = find(r_d_vec<= rs,1,'first');
[~,mui]         = find(mu_vec>= mu_ss_target,1,'first');
i_m_ss          = rs - 1/2*(1-varrho)*((1+omega)*chi_p_vec(mui)+(1-omega)*chi_m_vec(mui));
r_m_aux         = i_m_ss; 

% Interbank Market Equilibrium
MPCC_bankblock_ii;
MPCC_interbank_vecs;


zlb_index     = find((r_d_vec<=0),1,'first');
r_d_vec_aux   = r_d_vec;
r_b_vec_aux   = r_b_vec;
Dr_zlb        = Dr(zlb_index);

if zlb_index<length(mu_vec);
    r_d_vec_aux(zlb_index:end)=r_d_vec_aux(zlb_index);
    r_b_vec_aux(zlb_index:end)=r_b_vec_aux(zlb_index);
    Dr(zlb_index:end)=Dr(zlb_index);
end

figure(cc)
plot(mu_vec(mu_index),[r_b_vec_aux(mu_index)' r_d_vec_aux(mu_index)'],'LineWidth',3); hold on;
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

figure(cc)
plot(mu_vec(mu_index),Dr(mu_index)','LineWidth',3); hold on;
scatter(mu_ss_target, rsp_ss ,40,'MarkerEdgeColor','r','MarkerFaceColor','r');
label_y('$\Delta r_t$'); label_x('Liquidity Ratio $\Lambda$'); grid on; hold on;
axis tight;
formataxis(gca);
if printit==1
    imprime(['fig' nameplot num2str(cc)]);
    imprpdf(['fig' nameplot num2str(cc)]);
end
cc=cc+1;

% Interbank market plots
if plotib==1
    MPCC_interbank_plots;
end


