function [DV]=objfcn_multiple_SK(var,mu_tbp,mu_v,R_v,J2_v,x_0,t_0,N_orbits)
% INPUTS
% var - [14*N_orbitsx1] design variables
% x_0 - [6x1] initial state at the pericentre (fixed)
% N_orbits [1x1] number of orbits to perform the sk
% OUTPUTS
% DV - overall DV of the transfer
%
% Authors:
% Fabrizio Maccari, Matteo Lusvarghi
% 05-01-23 - first implementation
% -------------------------------------------------------------------------

options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);
%algorithm
%1st propagation
x_sk1=var(1:6);
t_sk1=var(7);
[~,prop_state] = ode113(@SCR3BP_dyn,[t_0 t_sk1],x_0,options_ode,mu_tbp,mu_v,R_v,J2_v);
v_prop = prop_state(end,4:6)';

%first DV contribution
DV=norm(v_prop-x_sk1(4:6));

for k=1:2*N_orbits-1
    %propgation
    t1=var(7*k);
    t2=var(7*(k+1));
    [~,prop_state] = ode113(@SCR3BP_dyn,[t1 t2],var(7*k-6:7*k-1),options_ode,mu_tbp,mu_v,R_v,J2_v);
    v_prop = prop_state(end,4:6)';
    %DV computation
    DV=DV+norm(v_prop-var(7*(k+1)-3:7*(k+1)-1));
end

end


