function [fb_ToF, timegrid, xstate, Sb, SAA] = EvalFlyBy(vinf_SOI, rSOI, rp, mu, X_p)
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% Date, User, brief summary of the modification
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
%

%% Function code

% Allocate input
v_inf_minus = vinf_SOI(1:3);
v_inf_plus = vinf_SOI(4:6);
R_p = X_p(1:3);
% V_p = X_p(4:6);


% Find Normal to the plane of the hypebola
u = cross(v_inf_minus,v_inf_plus); 
H_dir = u/norm(u);

% Find turning angle
% delta = (acos(dot(v_inf_minus/norm(v_inf_minus), v_inf_plus/norm(v_inf_plus))));

% Find velocity (scalar) at pericentre
v_p_min = sqrt(norm(v_inf_minus)^2 + 2*mu/rp);
v_p_plus = sqrt(norm(v_inf_plus)^2 + 2*mu/rp);

if abs(v_p_min - v_p_plus) > 1e-12
    warning("Pericentre velocities do not match (tol = 1e-8). Error = " + num2str(abs(v_p_min - v_p_plus)) + " km/s")
end

% Eccentricity of legs
% e1 = 1 + rp*norm(v_inf_minus)^2/mu_Planet;
% e2 = 1 + rp*norm(v_inf_plus)^2/mu_Planet;

% Semi major axis
% a1 = rp/(1-e1);
% a2 = rp/(1-e2);

% Angular momentum
% H1 = rp * v_p_min * [0; 0; 1];
% H2 = rp * v_p_plus * [0; 0; 1];
% beta = (pi - delta)/ 2;

% r_SOI = astroConstants(23) * 145.3; %% NB only for Earth (mass of the planet needed)

% fun1 = @(th) R_SOI - norm(H1)^2./mu * 1./(1+e1*cos(th));
% option = optimset('Display', 'off', 'TolX', 1e-6);
% th_inf = acos(-1/e1);
% 
% theta1 = fzero(fun1, 1.642, option);
% 
% fun2 = @(th) r_SOI - norm(H2)^2/mu * 1/(1 + e2 *cos(th));
% th_inf = acos(-1/e2);
% 
% theta2 = fzero(fun2, 1.642, option);
% 
% DT1 = evaluate_time_on_Hyperbola(theta1, -a1, mu_Planet, e1);
% DT2 = evaluate_time_on_Hyperbola(theta2, -a2, mu_Planet, e2);

% ToF_flyby = DT1 + DT2;

% Determine Initial condition for Incoming Hyperbola in Planet RF
DV_planet = v_inf_minus - v_inf_plus;
Rp_dir = DV_planet./norm(DV_planet);

RI = rp * Rp_dir;

VI = v_p_min * cross(H_dir, Rp_dir);
x0 = [RI'; VI']; % Initial conditions to propagate for incoming hyerbola

% Rotation_matrix1 = [cos(beta)+u(1)^2*(1-cos(beta)), u(1)*u(2)*(1-cos(beta))-u(3)*sin(beta), u(1)*u(3)*(1-cos(beta))+u(2)*sin(beta);
%     u(2)*u(1)*(1-cos(beta))+u(3)*sin(beta), cos(beta)+u(2)^2*(1-cos(beta)), u(2)*u(3)*(1-cos(beta))-u(1)*sin(beta);
%     u(1)*u(3)*(1-cos(beta))-u(2)*sin(beta), u(3)*u(2)*(1-cos(beta))+u(1)*sin(beta), cos(beta)+u(3)^2*(1-cos(beta))];

% RI = Rotation_matrix1 * RI; % Rotation radius of perigee from x-axis to injection ring 
% VI = cross(H1, RI); % Velocity direction
% VI = VI/norm(VI) * v_p_min; % Velocity vector at pericentre of incoming hyperbola


% VI = cross(H2,RI);
% VI = VI/norm(VI)*v_p_plus; % Velocity vector at pericentre of outgoing hyperbola
% y0_Hyp2 = [RI;VI];% Initial conditions to propagate for outgoing hyerbola
opts = odeset('RelTol', 3e-14, 'AbsTol', 3e-14, 'Events', @(t, x) SOIexit(t, x, rSOI));

dt = 60; % [s]
[timegrid, xstate, exit_time, ~] = ode113(@(t, x) RHS_2BP(t, x, mu), 0:dt:365.5*24*3600, x0, opts);
[~, xstate_back, ~, ~] = ode113(@(t, x) RHS_2BP(t, x, mu), -(0:dt:365.5*24*3600), x0, opts);

xstate = [flip(xstate_back(2:end, :)); xstate];

fb_ToF = 2*exit_time/3600; % [hours]

% Define r-theta-h frame at each epoch
r_state = xstate(:, 1:3);
% v_state = xstate(:, 4:6);

r_dir = r_state./vecnorm(r_state, 2, 2);
% v_dir = v_state./vecnorm(v_state, 2, 2);

% H_dir = cross(r_dir, v_dir);
H_dir = H_dir.*ones(length(r_dir), 3);
theta_dir = cross(H_dir, r_dir);

% Define SC position in SOI wrt SUN
R_SC = R_p' + r_state;

% Define Sb vector by projection at each epoch
S = R_SC./vecnorm(R_SC, 2, 2); % Position vector of the planet assumed equal to Sun rays incoming direction

Sb_r = dot(S, r_dir, 2);
Sb_theta = dot(S, theta_dir, 2);
Sb_h = dot(S, H_dir, 2);

Sb = [Sb_r, Sb_theta, Sb_h];
% Compute Solar Aspect Angle in deg from component normal to trajectory plane
SAA = mean(asind(Sb_h)); 


%% Local functions
    function [value, isterminal, direction] = SOIexit(~, x, rSOI)

        direction = 0;
        isterminal = 1;
        % Evaluate event condition
        value = norm(x(1:3)) - rSOI;


    end

    function [dxdt] = RHS_2BP(~, x, mu)

        dxdt = [x(4:6);
            -mu/norm(x(1:3))^3 * x(1:3)];

    end
end
