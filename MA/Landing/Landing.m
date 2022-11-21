%% Landing problem: Collocation method
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(0,'defaultAxesFontSize', 16)

%%
clearvars; close all; clc 

% step1: vehicle par (TO CHANGE!!!)
Tmax = 4*110e-3;                %[kN] Maximum Thrust  !!!!! 4* %%%%%%%%%%%%%%%%
Isp = 228;                      %[s] Specific Impulse
g0 = 9.81*1e-3;                 %[km/s^2] acceleration constant
m0 = 95+75.4;                   %[kg] initial mass of lander (both Non Sampling Orbiter- Sampling Lander and S-S)
m_dry = 75.4;

% Enceladus par
Re = 251.1;                                             %[km] mean radius of Enceladus
mu = (6.67430e-11 *8.6e19 )*10^(-9) ;                   % [km^3/s^2] Enceladus gravitational constant


% Adimensionalization
DU = Re;                        %[km]
TU = 3600;                      %[s]
VU = DU/TU;                     %[km/s]
MM = m0;                        %[kg]
acc = DU/(TU)^2;                %[km/s^2]
FU = MM*acc;                    %[kN]
GM = DU^3/TU^2;                 %[km^3/s^2]

Tmax = Tmax/FU;
Isp = Isp/TU;
g0 = g0/acc;
m0 = m0/MM;
m_dry = m_dry/MM;
Re = Re/DU;
mu = mu/GM;

% organize par 
par(1) = Tmax;
par(2) = Isp;
par(3) = g0;
par(4) = mu;
par(5) = Re;

N = 100;

% INITIAL ORBIT : circular, polar (TO CHANGE!!)
h = 100/DU;
r_mod = Re+h;
xi = -r_mod/sqrt(2);
yi = -r_mod/sqrt(2);
v_mod = sqrt(mu/norm([xi; yi]));
vxi = v_mod/sqrt(2);
vyi = -v_mod/sqrt(2);
state_i = [xi yi vxi vyi m0];

% NLP vars (x1, u1, ..., xN, uN, t1, tN)
step_st = length(state_i);           % 5: rr,vv,m
step_var = 5+3;                      % 8: rr,vv,m,u,ax,ay

%shift to satisfy boundaries
eps = 1e-6;

% t1 = 0/TU;                              %initial time guess
t1 = eps;
tN = 0.6*3600/TU;                          %final time guess
h = (tN-t1)/(N-1);
s0 = state_i;
guess = zeros(step_var*N+2, 1);
guess(1:5,1) = s0;
guess(6) = 1-eps;
guess(7:8) = -s0(3:4)'/norm(s0(3:4));
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

tspan_l(1) = t1;
tspan_l(N) = tN;
for k = 2:N-1
    tspan_l(k) = t1 + h*(k-1);
end

uk = guess(6:8);
for k = 1:N-1
    [~, output] = ode113(@landing_dyn, [tspan_l(k) tspan_l(k+1)], s0', options, uk, par);
    s0 = output(end,1:5);
    vv_vers = s0(3:4)'/norm(s0(3:4));
    if k == N-1
        uk = [1-eps; -vv_vers];
    else
        uk = [eps; -vv_vers];
    end

    guess(k*step_var+1:(k+1)*step_var) = [s0'; uk];
end

guess(end-1:end) = [t1; tN];

%one integration starting from second point
[~, output1] = ode113(@landing_dyn, [tspan_l(2) tspan_l(N)], guess(9:13)', options, [0;0;0], par);

% Check guess points
circ = 0:pi/1e4:2*pi;
r_i = state_i(1:2)';
v_i = state_i(3:4)';
t_range = [t1 15*tN];
[~, circ_orb] = ode113(@dyn, t_range, [r_i;v_i], options, par);
r_guess = circ_orb(:,1:2);
L = length(r_guess(:,1));
figure; hold on; grid on; grid minor
for k = 1:N
    plot(guess((k-1)*step_var+1), guess((k-1)*step_var+2), '.r', 'LineWidth', 1.2);
    plot(r_guess(:,1), r_guess(:,2), '-m', 'LineWidth', 0.5);
    plot(Re*exp(1i*circ),'-b')
    axis equal
end
plot(output1(:,1), output1(:,2), '.r', 'LineWidth', 1.2)
title('Initial Guess')
xlabel('$x$')
ylabel('$y$')

% check if initial guess satisfies constraints
for k = 1:N
    mass_check(k) = guess(step_var*(k-1)+5);
    u_check(k) = guess(step_var*(k-1)+6);
end
flag_mass = find(mass_check < m_dry);
flag_umin = find(u_check < 0);
flag_umax = find(u_check > 1);


% constraints for fmincon
A = [zeros(step_var*N,1); 1; -1]';
b = 0;
Aeq = [];
beq = [];

lb = -Inf*ones(1,step_var*N+2);
ub = Inf*ones(1,step_var*N+2);
for k = 1:N
    lb(step_var*(k-1)+5) = m_dry;
    lb(step_var*(k-1)+6) = 0;
    ub(step_var*(k-1)+6) = 1;
end
lb(end-1) = 0;
lb(end) = 0;

%%
% options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point',...
%     'SpecifyObjectiveGradient', true, 'MaxIter', 1000, 'MaxFunctionEvaluations', 1e5);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point',...
    'SpecifyObjectiveGradient', false, 'MaxIter', 1000, 'MaxFunctionEvaluations', 1e5);

[x_final, fval, exitflag, struct] = fmincon(@(var) land_objfun(var, state_i, par, N),guess,A,b,Aeq,beq,lb,ub, ...
    @(var) land_nonlincon(var, state_i, par, N),options);

%%
% Trajectory
figure; hold on; grid on; grid minor
for k = 1:N
    plot(x_final((k-1)*step_var+1), x_final((k-1)*step_var+2), '.r', 'LineWidth', 1.2);
    plot(r_guess(:,1), r_guess(:,2), '-m', 'LineWidth', 0.5);
    plot(Re*exp(1i*circ),'-b')
    axis equal
end
title('Optimized Landing Trajectory')
xlabel('$x$')
ylabel('$y$')


% Control law
figure; hold on; grid on; grid minor
t_plt = linspace(x_final(end-1), x_final(end), N);
for k = 1:N
    control(k) = x_final((k-1)*step_var+6);
    plot(t_plt(k),control(k), 'ob');
end
title('Control Law')
xlabel('$t$')
ylabel('$u$')

% Thrust
Thrust_min = min(control)*(Tmax*FU)*1e3;     %[N]

% Final Velocity
vv_fin = x_final(end-7:end-6)*VU;
v_fin = norm(vv_fin)*1e3;                    %[m/s]

% Propellant Mass
m_fin = x_final(end-5)*MM;
m0 = m0*MM;
m_prop = m0 - m_fin;

% Print Results
fprintf('RESULTS:\n\nMinimum thrust: %.4f [N]\nFinal velocity: %e [m/s]\nPropellant mass: %.1f [kg]', Thrust_min, v_fin, m_prop);


% TO DO:
% Change initial guess (?)
% Insert attitude

% ASK Guadagnini:
% Downrange as path constraint
% Initial guess: Do we need to divide landing if altitude is high?



