% Constructs Many Definition used in interbank market analysis
% For Interbank Market Spreads
fedbal_p   =@(mu) max(mu-varrho+omega*(1-varrho),0)                                   ;
fedbal_m   =@(mu) max(varrho-mu+omega*(1-varrho),0)                                   ; 
theta      =@(mu) fedbal_m(mu)./fedbal_p(mu)                                          ;
bartheta_go=@(theta) 1+(theta-1).*exp(barlam)                                         ;
bartheta_lo=@(theta) 1./bartheta_go(theta.^(-1))                                      ;
bartheta   =@(theta) (theta>=1).*bartheta_go(theta)+(theta<1).*bartheta_lo(theta)     ;

% Interbank market cost functions
% Definitions at interior
chi_d      =@(pol_d,theta) pol_d*(bartheta(theta)./theta).^(eta)                    ;
chi_p      =@(pol_d,theta) chi_d(pol_d,theta).*(theta.^(eta).*bartheta(theta).^(1-eta)-theta)./(bartheta(theta)-1);
chi_m      =@(pol_d,theta) chi_d(pol_d,theta).*(theta.^(eta).*bartheta(theta).^(1-eta)-1)./(bartheta(theta)-1);

% definitions at inf limit
chi_p_inf=@(pol_d) pol_d*(1-exp(-(1-eta)*barlam));
chi_m_inf=@(pol_d) pol_d                         ;
ibar_inf=@(pol_d) chi_p_inf(pol_d)/(1-exp(-barlam))              ;

% definitions at zero limit
chi_p_z=@(pol_d) 0                                                ;
chi_m_z=@(pol_d) pol_d*exp(-eta*barlam)                           ;
ibar_z =@(pol_d) pol_d*(1-(1-exp(-(eta)*barlam))/(1-exp(-barlam)));

% definitions at unit limit
chi_p_o=@(pol_d) pol_d*(1-eta)*(1-exp(-barlam));
chi_m_o=@(pol_d) pol_d*(1-eta*(1-exp(-barlam)));  
ibar_o =@(pol_d) pol_d*(1-eta)                 ;

% Definition of Rates
bal_p=@(mu) mu-varrho+(1-varrho)*omega;
bal_m=@(mu) mu-varrho-(1-varrho)*omega;
r_b_mu=@(r_m,mu,chip,chim) r_m+1/2*(chim.*((bal_p(mu)<0)+(bal_m(mu)<0))...
    +chip.*((bal_p(mu)>0)+(bal_m(mu)>0))); 
% At interior...
r_d_mu=@(r_m,mu,chip,chim) r_m+1/2*(1-varrho)*(chim*(1-omega)+(1+omega)*chip).*((bal_p(mu)>0).*(bal_m(mu)<0))+...
    (1-varrho)*(1+omega)*chip.*((bal_p(mu)>0).*(bal_m(mu)>0))+...
    (1-varrho)*(1-omega)*chim.*((bal_p(mu)<0).*(bal_m(mu)<0));

% Fed Profits per unit of Deposits - to get money value, multiply by M and
% mu inverse
dwprofs=@(pol_d,mu)     1/2*(chi_m(pol_d,theta(mu)).*fedbal_m(mu)-chi_p(pol_d,theta(mu)).*fedbal_p(mu));
dwprofs_z=@(pol_d,mu)   0;
dwprofs_o=@(pol_d,mu)   1/2*(chi_m_o(pol_d).*fedbal_m(mu)-chi_p_o(pol_d).*fedbal_p(mu));
dwprofs_inf=@(pol_d,mu) chi_m_inf(pol_d).*fedbal_m(mu);

% Fed Equity is expressed in terms of deposits
bsprofs=@(fedeq,r_b,r_m,mu) r_b.*fedeq+(r_b-r_m).*mu;

% Total profits given money multiplier and M
totprofs=@(bsprofs,dwprofs,M,mu) (dwprofs+bsprofs);

% Relevant Cutoffs
surplus_cut=(varrho-omega*(1-varrho))^(-1);
deficit_cut=(varrho+omega*(1-varrho))^(-1);