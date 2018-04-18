%% Construct Interest Rate Vector
r_vec=s_vec;
r_vec(s_vec>0)=rs_t(tt) ;
r_vec(s_vec<=0)=rb_t(tt);
