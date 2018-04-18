function [res]=MPCC_solve_steady_FP(rs,inputs,paths)
% Variables defined as global come in
MPCC_globals;

% Initial Guess for consumption vector
p=inputs.p;
V=inputs.V;

% Update Rate
% rs=exp(rs)                ; % allows only positive rates     
% rs=rs                       ; % impose borrowing limit

% [1] Solve for Steady State for each case
% Construct Interest Rate Vector
r_vec=s_vec       ; % Initialize it
r_vec(s_vec>0)=rs ;
r_vec(s_vec<=0)=rs+0.00; % pin down value

% Calls the Solver
MPCC_steady;

% Pick a clearing condition to solve
if strcmp(clearcond,'Y');
    res=Z_Y_ss;
elseif strcmp(clearcond,'S');
    res=Z_S_ss;
end