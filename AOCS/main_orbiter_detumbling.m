clearvars; close all; clc

DefaultFontSize = 16;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

rng shuffle

%% Model parameters
J_SO = diag([12043.36, 15990.19, 8547.84]); % [kg m^2]

muS = 3.7931187*1e16;            %[m^3/s^2]
muE = (6.6743e-11)*(1.0802e20);  %[m^3/s^2]

% Geometry
inc = deg2rad(45);  % Thruster inclination wrt face 
lx = 3.600; % [m]
ly = 2.500; % [m]
lz = 3.610; % [m]
D = 15*1e-2; % [m] Z offset wrt CM

J = J_SO;
ModelName = 'ADCS_orbiter';  


%% Initial conditions
wb0 = 0.001*rand(3, 1);
qrand = rand(4, 1);
qbn0 = qrand./norm(qrand);

% t0 = 0;
% tf = 1000;
% tspan = t0:0.1:tf;

% Science Orbit
n_days = 3;
filename = "Halo_" + num2str(n_days) + "days.mat";
state_str = load('Halo_3days.mat');
state = cell2mat(struct2cell(state_str))';
st_pos = state;
L = length(st_pos);

Xpos = reshape(st_pos, [6, 1, L]);

StopTime = n_days*24*3600;
period = 12*3600;

%% Sensors
% Gyros model parameters
% Nsamples = length(tspan);
sigma_RRW = 0.001;
sigma_ARW = 0.001;
% R_RRW = diag(ones(1, 3)*sigma_RRW^2);
% R_ARW = diag(ones(1, 3)*sigma_ARW^2);
% RRW_noise = mvnrnd(zeros(3, 1), R_RRW, Nsamples);
% ARW_noise = mvnrnd(zeros(3, 1), R_ARW, Nsamples);

%% Actuators
% 8 actuators 

% MR 103J 1N Rocket Engine Thrusters
% MIB = 0.0133 Ns
F_min = 0.19;  % N
F_max = 1.13; % N
MIB_real=F_min*0.1;

MIB_max=0.05;
Dt_min=MIB_max/F_min;

R = 0.5*[-ly*cos(inc)              ly*cos(inc)            ly*cos(inc)            -ly*cos(inc)             -ly*cos(inc)               ly*cos(inc)              ly*cos(inc)             -ly*cos(inc)          
          D*sin(inc)+lx*cos(inc)   D*sin(inc)-lx*cos(inc) D*sin(inc)+lx*cos(inc)  D*sin(inc)-lx*cos(inc)  -D*sin(inc)-lx*cos(inc)   -D*sin(inc)+lx*cos(inc)  -D*sin(inc)-lx*cos(inc)  -D*sin(inc)+lx*cos(inc)   
         -ly*sin(inc)             -ly*sin(inc)            ly*sin(inc)             ly*sin(inc)              ly*sin(inc)               ly*sin(inc)             -ly*sin(inc)             -ly*sin(inc)             ];

w = ones(length(R(1,:)),1);

R_ast = pinv(R);

%% Control 
k_p = 0;
k_d = 100;
w_c = [0 0 0]';

%% Run Simulink
% out = sim(ModelName, 'Timeout', StopTime);
return