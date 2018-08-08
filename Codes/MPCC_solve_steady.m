function [res]=MPCC_solve_steady(rs,inputs)
% Variables defined as global come in
MPCC_globals;

% Initial Guess for consumption vector
p=inputs.p;
V=inputs.V;

% [1] Solve for Steady State for each case
% Construct Interest Rate Vector
if any(strcmp(mpregime,{'FP' 'RP'}))
    r_vec=s_vec       ; % Initialize it
    r_vec(s_vec>0)=rs ;
    r_vec(s_vec<=0)=rs; % pin down value
elseif any(strcmp(mpregime,{'BP' 'MP'}))
    r_vec=s_vec       ; % Initialize it
    r_vec(s_vec>0)=rs ;
    r_vec(s_vec<=0)=rs+rsp_ss; % pin down value
end
% Calls the Solver
MPCC_steady;

% Pick a clearing condition to solve
if strcmp(clearcond,'Y');
    res=Z_Y_ss;
elseif strcmp(clearcond,'S');
    res=Z_S_ss;
end