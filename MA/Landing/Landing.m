%% Landing problem: Collocation method
clearvars; close all; clc

% step1: vehicle par (TO CHANGE!!!)
Tmax = 4*110e-3;                %[kN] Maximum Thrust  !!!!! 4* %%%%%%%%%%%%%%%%
Isp = 228;                      %[s] Specific Impulse
g0 = 9.81*1e-3;                 %[km/s^2] acceleration constant
m0 = 95+75.4;                   % [kg] initial mass of lander (both Non Sampling Orbiter- Sampling Lander and S-S)
m_dry = 75.4;

% Enceladus par
Re = 251.1;                                             %[km] mean radius of Enceladus
mu = (6.67430e-11 *8.6e19 )*10^(-9) ;                   % [km^3/s^2] Enceladus gravitational constant


% organize par 
par(1) = Tmax;
par(2) = Isp;
par(3) = g0;
par(4) = mu;
par(5) = Re;

N = 50; 

% INITIAL CONDITION : pericenter conditions of Enceladus Orbit
% xi = Re+400;   
% vyi = 0.1;
% INITIAL ORBIT : circular, polar, altitude = 20 km. (TO CHANGE!!)
xi = 0;
yi = Re+20;
vxi = -sqrt(mu/norm([xi; yi]));
vyi = 0;
state_i = [xi yi vxi vyi m0];

% NLP vars (x1, u1, ..., xN, uN, t1, tN)
step_st = length(state_i);           % 5: rr,vv,m
step_var = 5+3;                      % 8: rr,vv,m,u,ax,ay

t1 = 0;                              %initial time guess
tN = 30*60;                          %final time guess
h = (tN-t1)/(N-1);
s0 = state_i;
guess = zeros(step_var*N+2, 1);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
for k = 1:N
    tk = t1 + h*(k-1);
    t_next = tk + h;
    vv_vers = s0(3:4)'/norm(s0(3:4));
    if k == N || k == 1
        uk = [1; -vv_vers];
    else
        uk = [0; -vv_vers];
    end

    [~, output] = ode113(@landing_dyn, [tk t_next], s0', options, uk, par);
    % state_guess(k+1,1:5) = output(end,1:5);
    % s0 = state_guess(k+1,1:5);
    s0 = output(end,1:5);
    guess((k-1)*step_var+1:k*step_var) = [s0'; uk];
end
guess(end-1:end) = [t1; tN];

% Check guess points
r_i = state_i(1:2)';
v_i = state_i(3:4)';
t_range = [t1 7*tN];
[~, out_guess] = ode113(@dyn, t_range, [r_i;v_i], options, par);
r_guess = out_guess(:,1:2);
L = length(r_guess(:,1));
figure; hold on; grid on; grid minor
for k = 1:N
    plot3(guess((k-1)*step_var+1), guess((k-1)*step_var+2), 0, 'or', 'LineWidth', 1.2);
    plot3(r_guess(:,1), r_guess(:,2), zeros(L,1), '--r', 'LineWidth', 1.2);
    axis equal
end

return

% constraints for fmincon
A = [zeros(step_var*N,1); 1; -1]';
b = 0;
Aeq = [];
beq = [];

lb = -Inf*ones(1,step_var*N+2);
ub = Inf*ones(1,step_var*N+2);
for k = 1:N %%%%%%%%%%%%%%%%%%%%%%
    lb(step_var*(k-1)+5) = m_dry;
    lb(step_var*(k-1)+6) = 0;
    ub(step_var*(k-1)+6) = 1;
end
lb(end-1) = 0;  %%%%%%%%%%%%
lb(end) = 0;


%call solver
% options = optimoptions('fmincon', 'Display', 'iter', 'FunctionTolerance', 1e-10, ...
 %   'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'MaxIter', 1000, ...
%    'MaxFunctionEvaluations', 1000000);
options = optimoptions('fmincon', 'Display', 'iter', ...
   'MaxIter', 1000, ...
    'MaxFunctionEvaluations', 100000);
fmincon(@(var) land_objfun(var, state_i, par, N),guess,A,b,Aeq,beq,lb,ub, ...
    @(var) land_nonlincon(var, state_i, par, N),options);

% TO DO:
% plot final solution
% change FE with Hermite-Simpson and trapezoidal with Gauss
% change initial guess (?)
% increase N

