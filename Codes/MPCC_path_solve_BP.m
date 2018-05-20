function [res,aux,D_t]=MPCC_path_solve_BP(rs_t,init_con,steady,paths)
% Declare global variables
MPCC_globals;

% initial condition for distribution 
p_t=rs_t*0+1;
f_0=init_con.f_0;
rsp_t=paths.rsp_t;

% Initial Value function
V=steady.V_ss;
f_ss=steady.f_ss; % Used only to draw a picture for the convergence
D_ss  = steady.D_ss;
rs_ss = steady.rs_ss;
RT_ss = steady.RT_ss;
TE_ss= steady.TE_ss;

if ~isempty(Dr_zlb)
    rsp_t=max(rsp_t,Dr_zlb);
end

% Update Borrowing Rate
rb_t=rs_t+rsp_t;

% Backout Transfer set
T_t=paths.T_t                          ;
index_T=(T_t~=0)                       ;
index_T_T=find(index_T==1,1,'last')    ;
RT_t=RT_ss*ones(1,T)                   ;
RT_t_aux=RT_t                          ;
TE_t=TE_ss*ones(1,T)                   ;
TE_t_aux=TE_t                          ;

% Test Mixer
for tt=1:T
    if tt<=T_post
        RT_t(tt  ) = (1+T_t(tt))*RT_ss                  ;
        TE_t(tt+1) = (1+rs_t(tt))*TE_t(tt)-RT_t(tt)   ;
    elseif tt>T_post
        aux_weight=delta_trans^(tt-T_post)                    ;
        TE_t_aux(tt+1) = TE_ss+delta_ef*(TE_t(tt)-TE_ss)      ;
        RT_t_aux(tt) = (1+rs_t(tt))*TE_t(tt)-TE_t_aux(tt+1)   ; % Find transfer rule that clears things
        RT_t(tt)     = (1-aux_weight)*RT_t_aux(tt)...
                        +aux_weight*RT_t(T_post);
        TE_t(tt+1) = (1+rs_t(tt))*TE_t(tt)-RT_t(tt) ;
    end
end

% Solve again - with the new price
MPCC_path;

% Plot Iteration
Periods = (1:1:T);
if plotit==1
    if rand<0.01
        if rand<0.04
            close all;
        end
        figure(cc);
        subplot(3,1,1);
        plot(Periods, Z_M); 
        title('Excess Money Supply - (Convergence)'); ylabel('Residual over Level'); grid on;
        hold on; drawnow;
        subplot(3,1,2);
        plot(Periods, Z_Y); 
        title('Excess Goods Supply - (Convergence)'); ylabel('Residual over Level'); grid on;
        hold on; drawnow;
        subplot(3,1,3);
        plot(Periods, Z_S); 
        title('Excess Asset Supply - (Convergence)'); ylabel('Residual over Level'); grid on;
        hold on; drawnow;
        figure(cc+1);
        plot(rs_t); drawnow;
    end
end

% Compute Invariant distribution again
if strcmp(clearcond,'Y');
    res=Z_Y;
elseif strcmp(clearcond,'S');
    res=Z_S;
end
aux = Z_S.*D_t;

end
