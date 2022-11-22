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

load('NewGlobalMin22112022.mat')

% Define trajectory features (N of FBs and sequence)
N = 3; % NUMBER OF FBs

planets = {'Earth','Venus', 'Earth', 'Earth','Jupiter', 'Saturn'};

planets_id = [3,2,3,3, 5, 6];

% Target orbit at Saturn
Ra_target = 200*R_Saturn;
Rp_target = 3*R_Saturn;

%% Analyze solution

initial_guess = NLPoptset_local(33,:,3);

DV_opt = objfun_EarthSaturntransfer_plot(initial_guess, planets_id, planets, Ra_target, Rp_target);

%% Plot the space of the variables
n1 = length(squeeze(NLPoptset_local(:,1,1)));
n2 = length(squeeze(NLPoptset_local(1,:,1)));
n3_th = length(squeeze(NLPoptset_local(1,1,:)));
t1 = LB(1, 1);
% Find the actual number of iterations performed
n3 = n3_th;
for k = 1:n3_th
    if squeeze(NLPoptset_local(:,:,k)) == zeros(n1,n2)
        n3 = n3 - 1;
    end
end

ind_v = [1,N+2:(N+2+N),15:18, 19:22];
ind_names = {'Departure Time [years from t1]','$ToF_1$','$ToF_2$',...
    '$ToF_3$','$ToF_4$','ToF_5','rp_1','rp_2','rp_3','rp_4', 'beta1', 'beta2', 'beta3', 'beta4'};
count = 1;

for j = ind_v
    figure
    title(ind_names{count})
    ylabel('Cost [Km/s]')
    xlabel(ind_names{count})
    hold on
    grid minor
    for k = 1:n3
        variables = squeeze(NLPoptset_local(:,j,k));
        costs = squeeze(feval_local(:,1,k));
        if j == 1
            variables = (variables - t1)/365.5;
        end

        scatter(variables,costs,10,'filled','DisplayName',['Iter: ',num2str(k)])
    end
    legend()
    count = count + 1;
end

% Find total time of flight
total_tof = sum(NLPoptset_local(N+2:(N+2+N)))./365.5;
