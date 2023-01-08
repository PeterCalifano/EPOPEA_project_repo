function [c,ceq] = nlcon_EOL(var,mu_tbp,t_0,state0_Halo,R_Saturn,r_max,mu_v,R_v,J2_v)
% INPUTS
% var - [8 x 1] design variables
% R_Saturn [adimensional] - Adimensional Saturn Radii
% OUTPUTS
% c - non linear inequality constraints
% ceq non linear equality constraints
%
% Authors:
% Matteo Lusvarghi
% 07-01-23 - first implementation
% -------------------------------------------------------------------------

%c=zeros(3,1);
ceq=zeros(3,1);

options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);

% Unpack Variables
x1 = var(1:6);
t1 = var(7);
x2 = var(8:13);
t2 = var(14);
tf = var(15);

% Propagate the Halo until t1
[~,prop_state] = ode113(@SCR3BP_dyn,[t_0 t1],state0_Halo,options_ode,mu_tbp,mu_v,R_v,J2_v);
state_Halo_t1 = prop_state(end,:)';

ceq(1:3) = state_Halo_t1(1:3) - x1(1:3);

% Propagate the State until tf
[~,prop_state] = ode113(@SCR3BP_dyn,[t1 t2],x1,options_ode,mu_tbp,mu_v,R_v,J2_v);
x1_prop = prop_state(end,:)';

ceq(4:6) = x1_prop(1:3) - x2(1:3);

% Propagate the State until tf
[~,prop_state] = ode113(@SCR3BP_dyn,[t2 tf],x2,options_ode,mu_tbp,mu_v,R_v,J2_v);
xf = prop_state(end,:)';

c(1) = norm(xf(1:3) + [mu_tbp;0;0]) - R_Saturn;
c(2) = norm(x2(1:3) + [mu_tbp;0;0]) - r_max;



end

