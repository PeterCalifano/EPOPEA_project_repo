function [dxdt] = landing_dyn (~, var, par)
% UNFOLD:
% state: r, v, m
r = var(1:2);
v = var(3:4);
m = var(7);

% par: Tmax, Is, g0
Tmax = par.Tmax;
Isp = par.Isp;
g0 = par.g0;
mu = par.mu;

% control u 
u = var(8);
alpha = var(9:10);

% dynamics 
dxdt = [v; 
    -mu/(norm(r)^3).*r + (Tmax/m)*u*alpha
    -Tmax*u/g0/Isp];

end