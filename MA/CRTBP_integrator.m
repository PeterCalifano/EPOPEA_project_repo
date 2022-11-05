function [output] = CRTBP_integrator(~, state, mu, dfdt)
% 2D Restricted Three Body Problem Integrator
% INPUT
%   state [42x1] first 6 elements represent the state while the others
%                represent te STM reshaped (from 6x6 to 36x1)
%   df     sym   symbolic matrix for the computation of A
%
% OUTPUT
%   
%------------------------------------------------------------------------------

% system state
x = state(1);
y = state(2);
z = state(3);
vx = state(4);
vy = state(5);
vz = state(6);

rr1 = [x+mu; y; z];
rr2 = [x+mu-1; y; z];
r1 = norm(rr1);
r2 = norm(rr2);

dxdt = [vx
        vy
        vz
        2*vy + x - (1-mu)*(x+mu)/(norm([x+mu; y; z]))^3 - mu*(x+mu-1)/(norm([x+mu-1; y; z]))^3
        -2*vx + y - (1-mu)*y/(norm([x+mu; y; z]))^3 - mu*y/(norm([x+mu-1; y; z]))^3
        -(1-mu)*z/(norm([x+mu; y; z]))^3 - mu*z/(norm([x+mu-1; y; z]))^3];

% STM
psi = state(7:end)';

% reshape psi
l = sqrt(length(psi));
psi_r = reshape(psi, l, l);

A = dfdt(x,y,z);

%A = double(subs(dfdt));
psi_dot = A * psi_r;
dpsi = reshape(psi_dot, 1, l*l)';

output = [dxdt; dpsi];

end