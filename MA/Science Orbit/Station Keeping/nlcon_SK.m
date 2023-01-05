function [c,c_eq,Gc, Gc_eq] = nlcon_SK(var,mu_tbp,mu_v,R_v,J2_v,r_nom_i,r_nom_f,threshold_r)
% 
% Description:
%   The function implements the non-linear constraints for a 2D Simple
%   Shooting problem 
%
% Inputs:
% var = [x,y,vx,vy,t0,tf]- Variables vector
% r_nom_f [3 x 1] - Adimensional position in the 3BP frame at the next
%       control point.
% threshold_r - Upper bound of the adimensional distance between the
%           propagated state and the reference one
%
% Outputs:
% c - Vector of inequality constraints
% c_eq - Vector of equality constraints
% Gc - Gradient of inequality constraints
% Gc_eq - Gradient of equality constraints
%
% Author:
% MATTEO LUSVARGHI
% -------------------------------------------------------------------------

% Extract the variables
xx_0 = var(1:6);
t0 = var(7);
tf = var(8);

% Set inequality constraints to zero
c_eq = [];
Gc = [];
Gc_eq = [];

% Propagate 
options_ode = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13);
[~,prop_state] = ode113(@SCR3BP_dyn,[t0, tf],xx_0,options_ode,mu_tbp,mu_v,R_v,J2_v);
xx_f = prop_state(end,:)';

% Compute the inequality constraints
c(1) = norm(xx_0(1:3) - r_nom_i) - threshold_r;
c(2) = norm(xx_f(1:3) - r_nom_f) - threshold_r;

end