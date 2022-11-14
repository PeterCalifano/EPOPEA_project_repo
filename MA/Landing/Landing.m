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
par.Tmax = Tmax;
par.Isp= Isp;
par.g0 = g0;
par.m0 = m0;
par.mu= mu;
par.Re = Re;
par.m_dry = m_dry;


% INITIAL CONDITION X0 = [x0 y0 xdot0 vdot0 m0]
% par.x0 = x0


% N = 120; !!!MIT MSc THESIS

% step: initial conditions (xi, ui, ti, tf).

% lower boundary: m(t) >= m_dry


% step: inequality constraints --> Re^2 <= r^2

%% Function

% 2. Objective function


function  [F, JF] = ObjFun(variables,~)


t_1 = 0; % ???????
t_N = variables(end);
% Unpack the control

N = floor(length(variables)/8);
F = 0;
h = (t_N-t_1)/(N-1);
for k = 1:(N-1)
   u_k = variables(6+8*(k-1));
   u_next = variables(6+8*(k));
   F = F+(u_k+u_next)*h/2;                                                 % trapezoidal. change in Gauss!
end
JF = [];
% add derivative
%     if nargout>1
%        
% 
% 
%     end

end

% 3. Constraints

function [ C, Ceq, JC, JCeq ] = Constr(variables, par)
t_1 = 0; % ???????
t_N = variables(end);
N = floor(length(variables)/8);
h = (t_N-t_1)/(N-1);


% Retrieve par

Re = par.Re;
x0 = par.x0;
m_dry = par.m_dry;


Ceq = zeros(6*N+2,1);

rr_end = variables(end-7:end-6);
r_end = norm(rr_end);
vv_end = variables(end-5:end-4);
v_end = norm(vv_end);
Ceq(end) = v_end;
Ceq(end-1) = r_end -Re;

for k = 1:(N-1)
   var =  variables((k-1)*8+1:8*k);
   x_k = var(1:5);
   x_next = variables(k*8+1:k*8+5);
   t_k = t_1+h*k;
   f_k = landing_dyn(t_k, var, par);
   epsilon = x_next-x_k-h*f_k;
   Ceq(5*(k-1)+1:5*k) = epsilon;
   alpha_x = var(end-1);
   alpha_y = var(end);
   Ceq(5*N+k) = alpha_x^2 + alpha_y^2 -1;
end
Ceq(5*(N-1)+1:5*N) = variables(1:5)-x0; % initial condition (could be linear)

var =  variables((N-1)*8+1:8*N);
alpha_x = var(end-1);
alpha_y = var(end);
Ceq(4*N+N) = alpha_x^2 + alpha_y^2 -1;

% derivative
% if nargout > 2
% 
% 
% end

% NB : many could be linear inequalities

C = zeros(8*N+1,1);
C(end) = -t_N;
C(end-1) = t_1-t_N;

var = variables(1:8);
x = var(1); 
y = var(2);
u = var(6);
alpha_x = var(7);
alpha_y = var(8);
C(1) = Re^2-x^2-y^2;
C(2) = -u;
C(3) = -u-1;
C(4) = -alpha_x;
C(5) = -alpha_x-1;
C(6) = -alpha_y;
C(7) = -alpha_y-1;
for k = 2:N
    var = variables((k-1)*8+1:8*k);
    x = var(1); 
    y = var(2);
    m = var(5);
    u = var(6);
    alpha_x = var(7);
    alpha_y = var(8);
    C(8*(k-2)+8) = Re^2-x^2-y^2;
    C(8*(k-2)+9) = -u;
    C(8*(k-2)+10) = -u-1;
    C(8*(k-2)+11) = -alpha_x;
    C(8*(k-2)+12) = -alpha_x-1;
    C(8*(k-2)+13) = -alpha_y;
    C(8*(k-2)+14) = -alpha_y-1;
    C(8*(k-2)+15) = m_dry -m;
end
JC = [];
JCeq = [];
end