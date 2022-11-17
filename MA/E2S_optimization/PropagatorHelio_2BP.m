function [r, v, timevec] = PropagatorHelio_2BP(initial_state, tf,mu_SUN)
%% PROTOTYPE
% [time, r, v] = twobp_solver(T, initial_state, mu, perturbation, dT)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Solves the two body problem for a given planetary constant mu,
% initial position r and velocity v, using ode113
% -------------------------------------------------------------------------
%% INPUT
% initial_state: [6x1] Initial state of the orbiting body (rx ry rz vx vy vz) [L] [L/T^2]
% -------------------------------------------------------------------------
%% OUTPUT
% r: [Nx3] position of the body [L] 
% v: [Nx3] velocity of the body [L/T^2]
% time: [Nx1] vector of time instants
% Note: N is the number of points representing the orbit
% -------------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano
% -------------------------------------------------------------------------
%% CHANGELOG
% v1: 2BP solver for Heliocentric trajectories, 17/11/2022
% -------------------------------------------------------------------------
%% Next upgrades
%mu_SUN = 1.32712440017987e+11;

%% Function code

r_init = initial_state(1:3);
v_init = initial_state(4:6);


options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
sysF = @(t, state_vec) [state_vec(4:6);...
    (-mu_SUN/norm(state_vec(1:3)).^3)*state_vec(1:3)];

[timevec, output_state] = ode113(sysF, [0,tf], [r_init; v_init], options);

r = output_state(:, 1:3)';
v = output_state(:, 4:6)';

end