function [chi_p_vec,chi_m_vec,btheta_vec,i_bar,psi_p,psi_m]=interbank_market(techpars,polpars,theta_vec);

% Takes values from 0 to inf in theta, technology parameters and policy
% parameters
eta=techpars.eta      ; % bargaining parameter
barlam=techpars.barlam; % efficiency parameter
pol_d=polpars.pol_s   ; % policy spread

% First Compute theta at limit
bartheta_go=@(theta) 1+(theta-1).*exp(barlam)                                       ;
bartheta_lo=@(theta) 1./bartheta_go(theta.^(-1))                                    ;
bartheta   =@(theta) (theta>=1).*bartheta_go(theta)+(theta<1).*bartheta_lo(theta)   ;

% Interbank Market Rates
chi_d      =@(pol_d,theta,btheta) pol_d*(btheta./theta).^(eta)                      ;

% Analytic Solutions
chi_p_vec=1:length(theta_vec)*NaN;
chi_m_vec=1:length(theta_vec)*NaN;
for th=1:length(theta_vec);
    theta=theta_vec(th)   ;
    btheta=bartheta(theta);
    btheta_vec(th)=btheta ;
    
    % Breaking the formula
    if theta==inf 
        chi_p      =pol_d*(1-exp(-barlam))  ;
        chi_m      =pol_d                   ;
    elseif theta==1
        chi_p      =(1-eta)*exp(-barlam)    ;
        chi_m      =eta*exp(-barlam)        ;
    elseif theta==0
        chi_p      =0;
        chi_m      =pol_d*eta*exp(-barlam);
    else
        chi_p=-chi_d(pol_d,theta,btheta).*(((btheta./theta).^(-eta)-1).*theta)./(btheta-1);
        chi_m=chi_d(pol_d,theta,btheta).*((btheta./theta).^(-eta).*theta-1)./(btheta-1)  ;
    end
    chi_p_vec(th)=chi_p;
    chi_m_vec(th)=chi_m;
end
% if theta=1 use