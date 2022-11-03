% Script prova
clear; clc; close all

% From ROM + margins due to stochastic
m0 = 1309*2;
dV = 3000; % trying Orbilander [m/s]

Is = 364;       % MMH, NTO
mp = preliminary_prop_mass(dV,m0,Is);

% MMH and NTO
propellant.rho_fu = 880;                 % [kg/m3] see reference on word
propellant.rho_ox = 1443;                % [kg/m3]
propellant.ratio = 2.50; 
% tanks
prop_tank.sigma = 950*1e+6; prop_tank.rho = 2780;     % from PS slides
gas_tank.sigma = 950*1e+6; gas_tank.rho = 2780;
% pressurant 
pressurant.R = 2077.3;                  % [J/kgK]
pressurant.gamma = 1.667;
pressurant.B = 3.5;

pressurant.mode = 'Blowdown';
[propellant_tank,pressurant_tank,m_gas] = biprop_sizing(mp,propellant,prop_tank,pressurant,gas_tank);