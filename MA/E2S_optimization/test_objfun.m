close all;
clear;
clc;

TU = 365*30*24*3600;
DU = astroConstants(1);
R_Saturn = astroConstants(26);

% Fix departure time 
t1 = 475848/TU;

% Initial guess vector 
var = [t1, t1+30*3600*24/TU, t1 + 70*3600*24/TU, DU, pi/8, 1000, t1 + 10*3600*24/TU,...
     DU, pi/4, 900, t1 + 50*3600*24/TU];

% Test of the objective function
DV = objfun_EarthSaturntransfer(var,1,[3,2,3],20*R_Saturn,200*R_Saturn,TU);

% Test with solver
stoptime = 100;
opts_ga = optimoptions('ga', 'FunctionTolerance', 1e-12, 'MaxTime', stoptime,...
    'UseParallel', true, 'PopulationSize', 1000, 'Display', 'iter', 'MaxGenerations', 1e6);

[opt_state, feval, exitflag, output] = ga(@(var) objfun_EarthSaturntransfer(var, 1, [3, 2, 3], 20*R_Saturn,200*R_Saturn,TU),...
    length(var), [], [], [], [], [], [], [], [], opts_ga);
