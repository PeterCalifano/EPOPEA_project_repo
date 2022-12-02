close all
clear
clc

rng shuffle
%% Model parameters
J_NSO = diag([13760.61, 15399.52, 15621.60]); % [kg m^2]
J_SO = diag([16903.03, 17415.25, 24963.03]); % [kg m^2]

model = 1;
% Map: 1) NSO, 2) SO, 3) SL

switch model
    case 1
        J = J_NSO;
        ModelName = 'ADCS_orbiter';
    case 2
        J = J_SO;
        ModelName = 'ADCS_orbiter';
    case 3
       
end

%% Initial conditions
wb0 = 0.01*rand(3, 1);
qrand = rand(4, 1);
qbn0 = qrand./norm(qrand);

t0 = 0;
tf = 1000;
tspan = t0:0.1:tf;

% Gyros model parameters

% Random Noise Generation - Gyros
Nsamples = length(tspan);

sigma_RRW = 0.001;
sigma_ARW = 0.001;

R_RRW = diag(ones(1, 3)*sigma_RRW^2);
R_ARW = diag(ones(1, 3)*sigma_ARW^2);

RRW_noise = mvnrnd(zeros(3, 1), R_RRW, Nsamples);
ARW_noise = mvnrnd(zeros(3, 1), R_ARW, Nsamples);

% Misalignment
S = eye(3);

out = sim(ModelName, 'Timeout', tf);

