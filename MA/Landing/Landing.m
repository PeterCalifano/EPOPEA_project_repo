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
par.mu= mu;
par.Re = Re;

N = 120; 

% INITIAL CONDITION x0 = [x0 y0 vx0 vy0 m0]
state_0 = [x0 y0 vx0 vy0 m0];

% NLP vars (x1, u1, ..., xN, uN, t1, tN)


% lower boundary: m(t) >= m_dry
% m_dry = par.m_dry;

% step: inequality constraints --> Re^2 <= r^2

%% Function

% 2. Objective function




% 3. Constraints

