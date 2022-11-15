%% MISCELLANEOUS/TEST CODES - EARTH-SATURN OPTIMIZATION
clear; close all; clc

%% Load SPICE Kernels
cspice_kclear();
cspice_furnsh('..\..\spice_kernels/pck00010.tpc')
cspice_furnsh('..\..\spice_kernels/naif0012.tls')
cspice_furnsh('..\..\spice_kernels/gm_de431.tpc')
cspice_furnsh('..\..\spice_kernels/de440s.bsp')
%% DATA:
% Radius of Saturn from Spice
R_Saturn = astroConstants(26); % [km]

% Normalization coefficients
TU = 365*24*3600; % 1 year
DU = 1.495978707e+8; % 1 AU

% Initial Time manipulation
date_1 = '2030-01-01 00:00:00.00 UTC'; % First available date to launch
t1 = cspice_str2et(date_1)/TU;
t1 = t1/TU;

% Define trajectory features (N of FBs and sequence)
N = 3; % NUMBER OF FBs
planets = {'Earth','Venus','Earth','Jupiter Barycenter','Saturn Barycenter'};
planets_id = [3,2,3,5,6];

% Target orbit at Saturn
Ra_target = 200*R_Saturn;
Rp_target = 20*R_Saturn;

%% OPTIMIZATION
clc

% Initial guess

% Define linear constraints

% - t_1 < t_initial = (2030)








A = [];
b = [];

% Optimizer
DV = fmincon(@(var) objfun_EarthSaturntransfer(var,N,planets,Ra_target,Rp_target,TU),guess_0,A,b);



