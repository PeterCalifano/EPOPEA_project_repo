function [c,ceq] = nlcon_EOL(var,mu_tbp,x0_Halo)
% INPUTS
% var - [8 x 1] design variables
% 
% OUTPUTS
% c - non linear inequality constraints
% ceq non linear equality constraints
%
% Authors:
% Matteo Lusvarghi
% 07-01-23 - first implementation
% -------------------------------------------------------------------------

options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);

c=zeros(4*N_orbits,1);
ceq=zeros(6*N_orbits,1);

options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);

% Unpack Variables
x1=var(1:6);
t1=var(7);

% Propagate the Halo until t1
[~,prop_state] = ode113(@CR3BP_dyn,[t_0 t1],x0_Halo,options_ode,mu_tbp);
state_Halo_t1 = prop_state(end,:)';


end

