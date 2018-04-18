% Theta Vecs
fedbal_m_vec=fedbal_m(mu_vec)       ;
fedbal_p_vec=fedbal_p(mu_vec)       ;
theta_vec=theta(mu_vec)             ;
bartheta_vec=bartheta(theta(mu_vec));
chi_d_vec=chi_d(poldaux,theta_vec)  ;
chi_p_vec=chi_p(poldaux,theta_vec)  ;
chi_m_vec=chi_m(poldaux,theta_vec)  ;
dwprofs_vec=dwprofs(poldaux,mu_vec);

% Replacing Values at Extremes
chi_p_vec(theta_vec==inf)=chi_p_inf(poldaux);
chi_m_vec(theta_vec==inf)=chi_m_inf(poldaux);
chi_p_vec(theta_vec==1)=chi_p_o(poldaux);
chi_m_vec(theta_vec==1)=chi_m_o(poldaux);
chi_p_vec(theta_vec==0)=chi_p_z(poldaux);
chi_m_vec(theta_vec==0)=chi_m_z(poldaux);
dwprofs_vec(theta_vec==inf)=dwprofs_inf(poldaux,mu_vec(theta_vec==inf));
dwprofs_vec(theta_vec==1)=dwprofs_o(poldaux,mu_vec(theta_vec==1));
dwprofs_vec(theta_vec==0)=dwprofs_inf(poldaux,mu_vec(theta_vec==0));

% chi_m_vec(isnan(chi_m_vec))=poldaux*(exp(-barlam)-eta*(1-exp(-barlam)));
% chi_d_vec(isnan(chi_d_vec))=chi_m_vec(isnan(chi_d_vec))-chi_p_vec(isnan(chi_d_vec));
r_b_vec=r_b_mu(r_m_aux,mu_vec,chi_p_vec,chi_m_vec);
r_d_vec=r_d_mu(r_m_aux,mu_vec,chi_p_vec,chi_m_vec);
inf_index=(theta_vec==inf);
r_d_vec(inf_index)=mu_vec(inf_index)*r_m_aux+(1-mu_vec(inf_index)).*r_b_vec(inf_index)-dwprofs_vec(inf_index);

% Total Profits from FED
% Balance Sheet Profits
bsprofs_vec=bsprofs(fedeqaux,r_b_vec,r_m_aux,mu_vec);

% Total profits given money multiplier and M
totprofs_vec=totprofs(bsprofs_vec,dwprofs_vec,M_aux,mu_vec);

% Technology and Policy
techpars.eta=eta      ; % bargaining parameter
techpars.barlam=barlam; % efficiency parameter
polpars.pol_s=poldaux   ; % policy spread
[chi_p_vec_a,chi_m_vec_a,btheta_vec]=interbank_market(techpars,polpars,theta_vec);