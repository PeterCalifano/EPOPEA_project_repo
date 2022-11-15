%% Landing problem: Collocation method
clearvars; close all; clc

% step1: vehicle par (TO CHANGE!!!)
Tmax = 110e-3;                %[kN] Maximum Thrust
Isp = 228;                   %[s] Specific Impulse
g0 = 9.81*1e-3;             %[km/s^2] acceleration constant
m0 = 95+75.4;           % [kg] initial mass of lander (both Non Sampling Orbiter- Sampling Lander and S-S)
m_dry = 75.4;

% Enceladus par
Re = 251.1;                                             %[km] mean radius of Enceladus
mu = (6.67430e-11 *8.6e19 )*10^(-9) ;                   % [km^3/s^2] Enceladus gravitational constant


% organize par in a structure
par(1) = Tmax;
par(2) = Isp;
par(3) = g0;
par(4) = mu;
par(5) = Re;

N = 120; 

% INITIAL CONDITION : pericenter conditions of Enceladus Orbit
state_i = [xi yi vxi vyi m0];

% NLP vars (x1, u1, ..., xN, uN, t1, tN)
step_st = length(state_i);           % 5: rr,vv,m
step_var = 5+3;                      % 8: rr,vv,m,u,ax,ay

t1 = 0;                              %initial time guess
tN = 30*60;                          %final time guess
h = (tN-t1)/(N-1);
s0(1,5) = state_i;
for k = 1:N-1
    tk = t1 + h(k-1);
    t_next = tk + h;
    vv_vers = s0(3:4)/norm(s0(3:4));
    if k == N-1 || k == 1
        uk = [1, -vv_vers];
    else
        uk = [0, -vv_vers];
    end

    [~, output] = ode113(@landing_dyn, [tk t_next], s0, options, uk, par);
    s0(k,1:5) = output(end,1:5);

    
end


% constraints for fmincon
A = [zeros(step_var*N,1); 1; -1];
b = 0;
Aeq = [];
beq = [];

lb = -Inf*ones(1,step_var*N+2);
ub = Inf*ones(1,step_var*N+2);
for k = 1:length(lb)
    lb(step_var*(k-1)+5) = m_dry;
    lb(step_var*(k-1)+6) = 0;
    ub(step_var*(k-1)+6) = 1;
end
lb(end) = 0;

%call solver
options = optimoptions('fmincon', 'Display', 'iter', 'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'MaxIter', 1000, ...
    'MaxFunctionEvaluations', 2000);
fmincon(@(var) land_objfun(var, state_i, par, N),x0,A,b,Aeq,beq,lb,ub, ...
    @(var) land_nonlincon(var, state_i, par, N),options)

%% Function

% 2. Objective function




% 3. Constraints

