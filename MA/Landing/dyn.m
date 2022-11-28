function [dxdt] = dyn (~, state, par)
% UNFOLD
% state: r, v
r = state(1:3);
v = state(4:6);

% par: mu
mu = par(4);

% dynamics 
dxdt = [v; 
        -mu/(norm(r)^3).*r];

end