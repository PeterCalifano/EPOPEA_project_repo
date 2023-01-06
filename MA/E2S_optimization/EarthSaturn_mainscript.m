%% MISCELLANEOUS/TEST CODES - EARTH-SATURN OPTIMIZATION
clear; close all; clc

%% Load SPICE Kernels
cspice_kclear();
try
    cspice_furnsh('spice_kernels/pck00010.tpc')
    cspice_furnsh('spice_kernels/naif0012.tls')
    cspice_furnsh('spice_kernels/gm_de431.tpc')
    cspice_furnsh('spice_kernels/de440s.bsp')
catch
    cspice_furnsh('..\..\spice_kernels/pck00010.tpc')
    cspice_furnsh('..\..\spice_kernels/naif0012.tls')
    cspice_furnsh('..\..\spice_kernels/gm_de431.tpc')
    cspice_furnsh('..\..\spice_kernels/de440s.bsp')
end
%% DATA:
% Radius of Saturn from Spice
R_Saturn = astroConstants(26); % [km]
mu_S = astroConstants(4);

% Normalization coefficients
TU = 24*3600; % 1 day
DU = 1.495978707e+8; % 1 AU

% Initial Time manipulation
% date_1 = '2030-01-01 00:00:00.00 UTC'; % First available date to launch
% t1 = cspice_str2et(date_1);
% t1 = t1/TU;

% date_2 = '2050-12-31 00:00:00.00 UTC'; % First available date to launch
% t2 = cspice_str2et(date_2);
% t2 = t2/TU;

%load('NewGlobalMin22112022.mat')

%% Define trajectory features (N of FBs and sequence)
marker = 7;
% 1 --> E-VEJ-S
% 2 --> E-VEE-S
% 3 --> E-VEVE-S
% 4 --> E-VEEJ-S
% 5 --> E-VVE-S
% 6 --> E-J-S
% 7 --> E-EEJ-S
% 8 --> E-EJ-S

[planets_id,planets,N] = sequence_selector(marker);

% Target orbit at Saturn
Ra_target1 = 200*R_Saturn;
Rp_target1 = 2.55*R_Saturn;

%% Analyze solution
%load('E_EEJ_S_bestvalue1.016.mat');
load('E_EEJ_S_1.033_FBhigh.mat');
%load('E_EEJ_S_bestvalue1.14backup2035.mat');
%load('E_EEJ_S_bestvalue1.02backup2043.mat');

[m, index_iter] = min(min_at_iter(min_at_iter>0));
index_pos = min_pos(index_iter);

initial_guess = NLPoptset_local(index_pos, :, index_iter);

[DV_opt, DV_breakdown,T_FB] = objfun_EarthSaturntransfer_plot(initial_guess, planets_id, planets, Ra_target1, Rp_target1,'static');

T_FB = T_FB/(24*3600);
%%
dep_time = cspice_et2utc(initial_guess(1)*3600*24,'C',0 )
time1stdsm = cspice_et2utc((initial_guess(1) + initial_guess(9)*initial_guess(5))*3600*24,'C',0 )
time1stfb = cspice_et2utc((initial_guess(1) + sum(initial_guess(5)))*3600*24,'C',0 )
time2nddsm = cspice_et2utc((initial_guess(1) + sum(initial_guess(5)) + initial_guess(10)*initial_guess(6))*3600*24,'C',0 )
time2stfb = cspice_et2utc((initial_guess(1) + sum(initial_guess(5:6)))*3600*24,'C',0 )
time3rddsm = cspice_et2utc((initial_guess(1) + sum(initial_guess(5:6)) + initial_guess(11)*initial_guess(7))*3600*24,'C',0 )
time3stfb = cspice_et2utc((initial_guess(1) + sum(initial_guess(5:7)))*3600*24,'C',0 )
time4thdsm = cspice_et2utc((initial_guess(1) + sum(initial_guess(5:7)) + initial_guess(12)*initial_guess(8))*3600*24,'C',0 )
arr_time = 24*3600*(initial_guess(1) + sum(initial_guess(5:5+N)) );
arr_time = cspice_et2utc(arr_time,'C',0 )

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

N = length(planets_id) - 2;
ind_tof = 5:(5+N);
ind_rp = (2*N + 7):(2*N + 7 + N-1);
ind_beta = (3*N + 7) : (3*N + 7 + N - 1);

%%%%%% SELECT THE VARIABLES TO BE PLOTTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_v = [1,2,3,4];


ind_names_all = {'Departure Date','v inf','u','v','tof 1','tof 2','tof 3','tof 4',...
    'eta 1','eta 2','eta 3','eta 4','rp 1','rp 2','rp 3','beta 1','beta 2',...
    'beta 3','beta 4'};

ind_names = ind_names_all(ind_v);

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
            variables = (variables);% - t1)/365.5;
        end

        scatter(variables,costs,10,'filled','DisplayName',['Iter: ',num2str(k)])
    end
    legend()
    count = count + 1;
end

% Find total time of flight
total_tof = sum(NLPoptset_local(N+2:(N+2+N)))./365.5;

%%
for i = 1:10
    for k = 1:4
    cspice_et2utc(NLPoptset_local(i,1,k)*3600*24,'C',0 )
    end
end