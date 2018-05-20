% There is no so much difference yet.

y1_out = w1;
y2_out = w2;
cond = 2*tol;

while cond>tol
    % Call Solver HJB steady state solver
    y1=y1_out;
    y2=y2_out;
    y =[y1; y2];
    c_bl = [y1 + r_vec(1)*s_vec(1);y1*ones(N-1,1)];  
    
    MPCC_solve_cons_steady;

    ds=s_vec(2)-s_vec(1); 

    % [2] Call Solve Values targets
    % Solving for Steady-State Distributions
    [f_ss]=KFE_ss_implicit(muF_ss,muB_ss,sigma_ss,s_vec,Ib,If);

    % Aggregate endowment and consumption % - fix herehhhhhhhh
    Y_ss=sum(f_ss(:)*w2)*ds + sum(f_ss(s_bl_index)*ds)*(w1 - w2); % Careful with agents that are short
    C_ss=f_ss'*c_ss*ds                         ;

    % Check Goods clearing
    check_drifts=muF_ss'*f_ss*ds;

    % Savings+Borrowings=0
    S_ss=f_ss'*s_vec*ds                   ; % Nominal Money Stock

    % Compute Deposits
    index_o=find(s_vec>=0,1,'first');

    D_ss=f_ss(index_o:end)'*s_vec(index_o:end)*ds;
    M0_ss=varrho*D_ss;
    if index_o>1
        B_ss=-f_ss(1:index_o-1)'*s_vec(1:index_o-1)*ds;
    else
        B_ss=0;
    end
    
    % Update income process
    if strcmp(mpregime,'BP');
        RT_ss=rs*(Ef_ss*D_ss);
        FT_ss=rsp_ss*B_ss;
        y1_out=w1 + RT_ss + FT_ss;
        y2_out=w2 + RT_ss + FT_ss;
        y_out=[y1_out; y2_out];
        cond=max(abs(y-y_out)./y);
    elseif strcmp(mpregime,'RP');
        RT_ss=D_ss*(Ef_ss)*rs;
        y1_out=w1+RT_ss;
        y2_out=w2+RT_ss;
        y_out=[y1_out; y2_out];
        cond=max(abs(y-y_out)./y);
    elseif strcmp(mpregime,'FP');
        RT_ss=0 ;
        y1_out=w1 + RT_ss ;
        y2_out=w2 + RT_ss  ;
        y_out=[y1_out; y2_out];
        cond=max(abs(y-y_out)./y);
    end    
    
    
end

% Construct Clearing residuals
Z_Y_ss=(Y_ss-C_ss)/Y_ss;
if strcmp(mpregime,'FP');
    Z_S_ss=(B_ss-D_ss)/D_ss; % pin down value
elseif any(strcmp(mpregime,{'RP' 'BP'}));
    Z_S_ss=((1+Ef_ss)*D_ss-B_ss)/D_ss; % Market clearing condition        
end