%% MISCELLANEOUS/TEST CODES - EARTH-SATURN OPTIMIZATION
clear; close all; clc

%% Load SPICE Kernels
cspice_kclear();
cspice_furnsh('spice_kernels/pck00010.tpc')
cspice_furnsh('spice_kernels/naif0012.tls')
cspice_furnsh('spice_kernels/gm_de431.tpc')
cspice_furnsh('spice_kernels/de440s.bsp')
%% DATA:
% Radius of Saturn from Spice
R_Saturn = astroConstants(26); % [km]
mu_S = astroConstants(4);

% Normalization coefficients
TU = 24*3600; % 1 day
DU = 1.495978707e+8; % 1 AU

% Initial Time manipulation
date_1 = '2030-01-01 00:00:00.00 UTC'; % First available date to launch
t1 = cspice_str2et(date_1);
t1 = t1/TU;

date_2 = '2050-12-31 00:00:00.00 UTC'; % First available date to launch
t2 = cspice_str2et(date_2);
t2 = t2/TU;

% Define trajectory features (N of FBs and sequence)
N = 3; % NUMBER OF FBs

planets = {'Earth','Venus','Earth','Jupiter','Saturn'};

planets_id = [3,2,3,5,6];

% Target orbit at Saturn
Ra_target = 200*R_Saturn;
Rp_target = 3*R_Saturn;

%% OPTIMIZATION
clc

% Initial guess
guess_0 = [t1, 4, 0.2, 0.4,...
    200, 500, 3000, 5500,...
    0.4, 0.5, 0.6, 0.6,...
    1.7, 1.5, 1.7, ...
    pi/2, pi/9, pi/8];




tic
DV = objfun_EarthSaturntransfer_2(guess_0, planets_id, Ra_target, Rp_target);
toc
%% Define linear constraints

% - t_1 < t_initial = (2030)

A = [];
b = [];

% Optimizer
DV = fmincon(@(var) objfun_EarthSaturntransfer2(var,N,planets,Ra_target,Rp_target,TU),guess_0,A,b);

%% Analyze solution

initial_guess = NLPoptset_local(33,:,7);

DV_opt = objfun_EarthSaturntransfer_plot(initial_guess, planets_id, planets, Ra_target, Rp_target);








