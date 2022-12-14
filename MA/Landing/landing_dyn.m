function [dxdt] = landing_dyn (~, state, uvec, par)
% UNFOLD
% state: r, v, m
r = state(1:3);
v = state(4:6);
m = state(7);

% par: Tmax, Is, g0
Tmax = par(1);
Isp = par(2);
g0 = par(3);
mu = par(4);

% control u 
u = uvec(1);
alpha = uvec(2:4);

% dynamics 
dxdt = [v; 
        -mu/(norm(r)^3).*r + (Tmax/m)*u*alpha;
        -Tmax*u/(g0*Isp)];

end