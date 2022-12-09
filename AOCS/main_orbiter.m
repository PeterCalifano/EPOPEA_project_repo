close all
clear
clc

rng shuffle
%% Model parameters
J_NSO = diag([6277.21, 10241.83, 12845.78]); % [kg m^2]
J_SO = diag([6961.96, 16022.69, 18620.98]); % [kg m^2]

muS = 3.7931187*1e16;            %[m^3/s^2]
muE = (6.6743e-11)*(1.0802e20);  %[m^3/s^2]

% Map: 1) NSO, 2) SO, 3) SL
model = 1;

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

% Science Orbit
state_str = load('Halo_3days.mat');
state = cell2mat(struct2cell(state_str));
st_pos = state(:,:);
L = length(st_pos(:,1));
Xpos = reshape(st_pos', [6,1,L]);

period = 12*3600;

% % Gyros model parameters
% % Random Noise Generation - Gyros
% Nsamples = length(tspan);
% 
% sigma_RRW = 0.001;
% sigma_ARW = 0.001;
% 
% R_RRW = diag(ones(1, 3)*sigma_RRW^2);
% R_ARW = diag(ones(1, 3)*sigma_ARW^2);
% 
% RRW_noise = mvnrnd(zeros(3, 1), R_RRW, Nsamples);
% ARW_noise = mvnrnd(zeros(3, 1), R_ARW, Nsamples);
% 
% % Misalignment
% S = eye(3);


%out = sim(ModelName, 'Timeout', tf);

