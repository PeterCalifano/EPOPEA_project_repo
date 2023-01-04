function DV = objfun_SK(var,mu_tbp,mu_v,R_v,J2_v)

%%%%%%%%%%%%%%%%%%%%%%%%%%% TO BE COMPLETED %%%%%%%%%%%%%%%%%%%%%%%

% Unpack the variables vector
xx_0 = var(1:4);
t0 = var(5);
tf = var(6);

% Propagate 
options_ode = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13);
[~,prop_state] = ode113(@SCR3BP_dyn,[t0, tf],xx_0,options_ode,mu_tbp,mu_v,R_v,J2_v);
xx_f = prop_state(end,:)';


end