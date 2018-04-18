function [f_ss]=KFE_ss_implicit(muF,muB,sigma,s_vec,Ib,If)
% Computes Forward Kolmogorov Equation using implicit method

% Begin Construction of Matrices
N = length(s_vec); 
ds = s_vec(2:N) - s_vec(1:N-1);  

% Initializing volatility terms
S2 = zeros(N,1); 
S2(2:N-1) = sigma(2:N-1).^2./(ds(1:N-2) + ds(2:N-1));

% Off-Diagonal Elements
MU  = If.*muF;
DU = zeros(N,1); 
SU = zeros(N,1); % last
DU(1:N-1) = MU(1:N-1)./ds; % last
SU(1:N-1) = S2(2:N)./ds;     % last

% On Diagonal
MD  = -Ib.*muB;
DD      = zeros(N,1); 
SD = zeros(N,1);     % last
DD(2:N) = MD(2:N)./ds; % last
SD(2:N) = S2(1:N-1)./ds;     % last

% Diagonal Term
D0 = - (DU + DD + 2*S2/ds(1));

% Making it Compatible with Sparse Matrix Notation
DU=[0; DU(1:end-1)];
DD=[DD(2:end); 0];

SU=[0; SU(1:end-1)];
SD=[SD(2:end); 0];

% Building Matrix
A= spdiags([DD D0 DU],-1:1,N,N);
A=A';

B = spdiags([SD SU],[-1 1],N,N);

A = A+B;


b = zeros(N,1);
i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,N-i_fix)];
A(i_fix,:) = row;

f_ss = A\b;
f_sum = f_ss'*(ones(N,1).*[ds(1);ds]);
f_ss = f_ss./f_sum;

end