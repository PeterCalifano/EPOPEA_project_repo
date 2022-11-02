function [r, v, time] = twobp_solver(T, initial_state, mu, J2_effect, dT)
%% PROTOTYPE
% [time, r, v] = twobp_solver(T, initial_state, mu, perturbation, dT)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Solves the two body problem for a given planetary constant mu,
% initial position r and velocity v, using ode45
% "Default Scan rate" of the solver: 2 Hz
% -------------------------------------------------------------------------
%% INPUT
% T: [scalar] Final instant of time  [T] or "single rev"
% initial_state: [6x1] Initial state of the orbiting body (rx ry rz vx vy vz) [L] [L/T^2]
% mu: [scalar] Gravitational constant of the planet [L^3/T^2]
% J2_effect: [scalar] switch on (1) or off (0) the J2 effect
% dT: [scalar] Time step of the output
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
% v1: 2BP solver without perturbation, works correctly, coded before 30/10/2021
% v2: 2BP with J2 effect coded, 31/10/2021
% v3: modified T input to allow a single revolution propagation without
%     providing the period, 15/11/2021
% v4: backward propagation option implemented, 03/12/2021
% -------------------------------------------------------------------------
%% Next upgrades
% 1) Perturbed 2BP with additional perturbations
% 2) Check on relative magnenitudes of quantities and scale accordingly
% ONLY for execution

%% Function code
if ~exist('dT', 'var')
    dT = 0.5;
    disp(["Default time step: ", dT]);
end

if J2_effect == 0 || (nargin <= 3 && ~exist("perturbation", "var"))
    % Non perturbed solution
    disp("Non perturbed 2BP selected")
    r_init = initial_state(1:3);
    v_init = initial_state(4:6);

    if nargin == 3 || ischar(T) == 1
        a = -mu/(2*(0.5*dot(v_init, v_init) - mu/norm(r_init))); % from Mech energy equation
        T = 2*pi*sqrt(a^3/mu); % orbital period [1/s]
    end

    ic = [r_init; v_init];
    tspan = 0:dT:T; % "scan rate" of the solver: 2 Hz
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13 );
    sysF = @(t, state_vec) [state_vec(4:6);...
        (-mu/norm(state_vec(1:3))^3)*state_vec(1:3)];

    [time, output_state] = ode45(sysF, tspan, ic, options);

    r = output_state(:, 1:3);
    v = output_state(:, 4:6);

elseif J2_effect == 1

    % J2 effect (oblateness) included
    disp("Perturbed 2BP selected")
    r_init = initial_state(1:3);
    v_init = initial_state(4:6);

    if nargin == 3 % if T is not provided
        a = -mu/(2*(0.5*dot(v_init, v_init) - mu/norm(r_init))); % from Mech energy equation
        T = 2*pi*sqrt(a^3/mu); % orbital period [1/T]
    end

    R_E = astroConstants(23); % Km
    J2 = astroConstants(9); % Constant coefficient

    % r = nan(floor(T)+1, 3);
    % v = nan(floor(T)+1, 3);

    % Setting up the initial conditions to propagate
    r_0 = initial_state(1:3);
    v_0 = initial_state(4:6);

    % for i = 1:T+1

    %     r_norm = norm(r(i, :));

    %     aJ2 = 1.5 * ((J2*mu*R_E^2)./(r_norm^4)) * [(r(i, 1)/r_norm) * (5*(r(i, 3)/r_norm)^2 - 1);
    %                                                (r(i, 2)/r_norm) * (5*(r(i, 3)/r_norm)^2 - 1);
    %                                                (r(i, 3)/r_norm) * (5*(r(i, 3)/r_norm)^2 - 3)];

    ic = [r_0; v_0];

    tspan = 0:0.5:T;
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    [time, output_state] = ode45(@odesys_2BP_J2, tspan, ic, options);

    r = output_state(:, 1:3);
    v = output_state(:, 4:6);

end

end







