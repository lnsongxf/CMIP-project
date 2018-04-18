% Function that backs-up consumption given value function
% cF = Upinv(max(p*dvF,10^(-10)))           ; % First-order condition
% cB = Upinv(max(p*dvB,10^(-10)))           ; % First-order condition
cF = Upinv(p*dvF)           ; % First-order condition
cB = Upinv(p*dvB)           ; % First-order condition
% [cF(s_bl_index)]=min(cF(s_bl_index),c_bl(s_bl_index)); % include borrowing constraint limit
% [cB(s_bl_index)]=min(cB(s_bl_index),c_bl(s_bl_index)); % include borrowing constraint limit
