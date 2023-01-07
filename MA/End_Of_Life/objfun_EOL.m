function DV = objfun_EOL(var,mu_tbp,x0_Halo)
% INPUTS
% var - [8 x 1] design variables
% 
% OUTPUTS
% DV [adimensional]

options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);

% Unpack Variables
x1=var(1:6);
t1=var(7);

% Propagate the Halo until t1
[~,prop_state] = ode113(@CR3BP_dyn,[t_0 t1],x0_Halo,options_ode,mu_tbp);
state_Halo_t1 = prop_state(end,:)';

% Compute the EOL DV
DV = norm(state_Halo_t1(4:6) - x1(4:6));

end

