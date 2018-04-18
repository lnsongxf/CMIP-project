function [res]=MPCC_solve_steady_MP(rs,inputs,paths);
% Variables defined as global come in
MPCC_globals_MP;

% Initial Guess for consumption vector
% c=targets.c;
p=inputs.p;
c=inputs.c;

% Update Rate
% rs=exp(rs)                ; % allows only positive rates     
rs=rs                       ; % impose borrowing limit

% Construct Interest Rate Vector
r_vec=s_vec       ; % Initialize it
r_vec(s_vec>0)=rs ;
r_vec(s_vec<=0)=rb;

% Call Solver
MPCC_solve_steady;
% MPPC_value_cons_MP;

% Reporting Steady State Solutions
disp('**** Steady State Solution ***');
disp('Computation time:');
toc
disp(['Iterations:   ' num2str(iter)]);
toc

% [2] Call Solve Values targets
% Solving for Steady-State Distributions
[f_ss]=KFE_ss_implicit(mu_ss,sigma_ss,s_vec,dt);
f_ss=f_ss'                  ; 

% Aggregate endowment and consumption 
Y_ss=sum(f_ss(:)*w2) + f_ss(1)*(w1 - w2);
C_ss=f_ss'*c_ss                         ;

% Check Goods clearing

% Savings+Borrowings=0
S_ss=f_ss'*s_vec                   ; % Nominal Money Stock

% Compute Deposits
index_o=find(s_vec>=0,1,'first');

D_ss=f_ss(index_o:end)'*s_vec(index_o:end);
M0_ss=resreq*D_ss;
if index_o>1
    B_ss=-f_ss(1:index_o-1)'*s_vec(1:index_o-1);
else
    B_ss=0;
end

% Monetary Aggregates
M1_ss = D_ss        ;
MM_ss = M1_ss/M0_ss ;

if plotiter==1
    plot(s_vec,f_ss); xlim([s_vec(1) s_vec(end)]); ylim([0 max(f_ss)]); drawnow;
end
res=D_ss-B_ss;