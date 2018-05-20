%% MPCC_policypath
% script that updates variables used when solving a time path
rr   = paths.rr_t(tt);
w1   = paths.mu_w_t(1,tt); 
w2   = paths.mu_w_t(2,tt); 
s1   = paths.sigma_w_t(1,tt); 
s2   = paths.sigma_w_t(2,tt);
s_bl = paths.s_bl_t(tt);

% Update Rates
MPCC_update_rates;

% Update to include borrowing constraint
c_bl(1)    = y1 + r_vec(1)*s_vec(1); 
c_bl(2:N,1)= y1     ; % % check borrowing limit
