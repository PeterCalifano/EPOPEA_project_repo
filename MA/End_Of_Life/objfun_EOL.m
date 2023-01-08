function DV = objfun_EOL(var,mu_tbp,t_0,state0_Halo,mu_v,R_v,J2_v)
% INPUTS
% var - [8 x 1] design variables
% 
% OUTPUTS
% DV [adimensional]

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

% Compute the EOL DV
DV1 = norm(state_Halo_t1(4:6) - x1(4:6));

% Propagate the state from t1 to t2
[~,prop_state] = ode113(@SCR3BP_dyn,[t1, t2],x1,options_ode,mu_tbp,mu_v,R_v,J2_v);
state_biel = prop_state(end,:)';

DV2 =  norm(state_biel(4:6) - x2(4:6));

DV = DV1 + DV2;
end

