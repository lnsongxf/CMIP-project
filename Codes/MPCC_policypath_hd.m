%% MPCC_policypath
% script that updates variables used when solving a time path
r =paths.r_t(tt);
rs=paths.rs_t(tt); 
rb=paths.rb_t(tt); 
rd=paths.rd_t(tt); 
w1=paths.mu_w_t(1,tt); 
w2=paths.mu_w_t(2,tt); 
s1=paths.sigma_w_t(1,tt); 
s2=paths.sigma_w_t(2,tt);
s_bl=paths.s_bl_t(tt);

% Update to include borrowing constraint
c_bl(1)  = c_dl     ; 
c_bl(2:N,1)= w1     ; % % check borrowing limit

% Update Rates
MPCC_update_rates;