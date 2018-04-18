%% Construct Test of Interbank Market
% Policy and technology
varrho= 0.10 ; % Reserve Requirement
barlam= 1.5  ; % Efficiency interbank market
eta   = 0.2  ; % Bargaining
omega = 0.08 ; % Average size of shock - interbank

% Policy Parameters
imonmu_ss=varrho;
i_m_ss   =0.02  ;
isp_ss   =0.01  ;

% Vector of Variables
mu_vec=(0.0:0.01:1) 		;
i_m=0               	;
pol_s_vec=(0.0:0.01:1)      ;

% Three parameter model
techpars.eta=eta       	; % bargaining parameter
techpars.barlam=barlam 	; % efficiency parameter
polpars.pol_s=poldaux   ; % policy spread

% Definitions
MPCC_bankblock 			;

% Build Vector of Theta
theta_vec=theta(mu_vec) ;

% Intializing Matrices
chi_p_mat = NaN(length(pol_s_vec),lenght(mu_vec)) ;
chi_m_mat = NaN(length(pol_s_vec),lenght(mu_vec)) ;
btheta_mat= NaN(lenght(pol_s_vec),lenght(mu_vec)) ;
i_mat     = NaN(lenght(pol_s_vec),lenght(mu_vec)) ;
psi_p_mat = NaN(lenght(pol_s_vec),lenght(mu_vec)) ;
psi_m_mat = NaN(lenght(pol_s_vec),lenght(mu_vec)) ;
i_b_mat   = NaN(lenght(pol_s_vec),lenght(mu_vec)) ;
i_d_mat   = NaN(lenght(pol_s_vec),lenght(mu_vec)) ;

% First we construct a matrix keeping the corridor fixed, moving the interest on reserves
for ii=1:lenght(pol_s_vec)
		polpars.pol_s=pol_s_vec(ii)   ;
		[chi_p_mat(ii),chi_m_mat(ii),btheta_mat(ii),i_mat(ii),psi_p_mat(ii),psi_m_mat(ii)]=interbank_market(techpars,polpars,theta_vec);
		
end
mu_mat=ones(length(pol_s_vec),1)*mu_vec;
i_b_mat=r_b_mu(i_m,mu,chip_mat,chim_mat);
i_d_mat=r_b_mu(i_m,mu,chip_mat,chim_mat);



