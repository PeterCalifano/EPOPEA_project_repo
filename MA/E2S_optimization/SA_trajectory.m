close all
clear
clc

DefaultFontSize = 16;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

rng shuffle

%% Load SPICE
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

%% Load solution workspace
% Baseline
load('E_EEJ_S_bestvalue1.016.mat')

% Main backup
% load('E_EEJ_S_bestvalue1.14backup2035.mat')

% Secondary backup
% load('E_EEJ_S_bestvalue1.02backup2043.mat')

clearvars -except min_at_iter NLPoptset_local planet_seq min_pos

R_Saturn = astroConstants(26); % [km]
Ra_target = 200*R_Saturn;
Rp_target = 2.55*R_Saturn;

% Extract best solution
[~, index_iter] = min(min_at_iter(min_at_iter > 0));
index_pos = min_pos(index_iter);
Xsol = NLPoptset_local(index_pos, :, index_iter);

%% Build test distributions 
variance_selector = 1;

switch variance_selector
    case 0
        % Uncertainty coeff.
        coeff3 = 0.01; % Between 0 and 1
        % 1 Standard deviation defined as about 3% of the variable nominal value
        VarianceVec = (Xsol.*coeff./3).^2;
    case 1
        % Launch date [JD]
        sigma3_LD = 1; 

        % Vinf [km/s]
        sigma3_Vinf = 0.2;
        % Direction of Vinf [-]
        sigma3_u = 0.005;
        sigma3_v = 0.005;

        % Time of Flights of arcs [days]
        sigma3_ToF = [1, 1, 1, 1];

        % Fraction at which DSM occurs [0 to 1]
        sigma3_eta = [0.001, 0.001, 0.001, 0.001];

        % Pericentre radius in [Rplanets]
        sigma3_Rp = [0.0001, 0.0001, 0.00001];

        % Beta angle determining flyby plane [rad]
        sigma3_beta = deg2rad([2, 2, 2]);

        % Assembly of diagonal of Covariance
        VarianceVec = ([sigma3_LD, sigma3_Vinf, sigma3_u, sigma3_v, sigma3_ToF,...
            sigma3_eta, sigma3_Rp, sigma3_beta]/3).^2;
end

%%TO ADD: CHECK OF FEASIBILITY ON VARIABLES AFTER POOL GENERATION

% Define number of samples and population
Nsamples = 1e5;
Nvar = length(Xsol);
% Define matrix 
Sigma = diag(VarianceVec);
% Population pool samples randomly with Sigma covariance (independet r.v.)
Xpool = [Xsol; mvnrnd(Xsol, Sigma, Nsamples)];

%% Test population with direct evaluation of ObjFun
dataset = cell(Nsamples+1, 1);
DV = zeros(Nsamples+1, 1);

tic
parfor idp = 1:Nsamples+1
    try
        [DV(idp, 1), dataset{idp}] = objfun_EStransfer_analysis(Xpool(idp, :), planet_seq, Ra_target, Rp_target);
    catch
        DV(idp, 1) = nan;
        dataset{idp} = 'Error';
    end
    disp(['Processing sample n° ', num2str(idp)])
end

toc

figure('WindowState', 'maximized');
hold on;
scatter(1:Nsamples+1, DV, 10, 'o', 'LineWidth', 1, 'MarkerFaceColor', 'flat');
yline(2, 'k--', '2 km/s Line', 'FontSize', 14, 'LineWidth', 1.05);
grid minor;
axis auto;
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
xlabel('Sample id')
ylabel('DV Cost [km/s]')
ax.LineWidth = 1.08;

% Evaluate how many below thr = 2 km/s
Mask1 = DV < 2;
Mask2 = DV < 1.5;
Mask3 = DV < 1.2;
count = [sum(Mask1), sum(Mask2), sum(Mask3)];
fprintf('\nN° samples < 2 km/s: %2g/%2g\n', count(1), Nsamples+1);
fprintf('N° samples < 1.5 km/s: %2g/%2g\n', count(2), Nsamples+1);
fprintf('N° samples < 1.2 km/s: %2g/%2g\n', count(3), Nsamples +1);
fprintf('N° errors: %2g\n\n', sum(isnan(DV)));

% Plot Launch date vs DV
figure('WindowState', 'maximized');
hold on;
scatter(Xpool(:, 1), DV, 9, 'o', 'LineWidth', 1, 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 0.8);
scatter(Xsol(1), DV(1), 100, 'x', 'LineWidth', 1.6, 'Color', 'r');
yline(2, 'k--', '2 km/s', 'FontSize', 16, 'LineWidth', 1.1);
yline(1.5, '--', '1.5 km/s', 'Color', '#fc8803', 'FontSize', 16, 'LineWidth', 1.1);
yline(1.2, '--', '1.2 km/s','Color', '#44e366','FontSize', 16, 'LineWidth', 1.1);
yline(1, '--', '1 km/s', 'Color', '#4444bb', 'FontSize', 16, 'LineWidth', 1.1);
grid minor;
axis auto;
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
title('Launch date vs DV cost')
xlabel('Launch date')
ylabel('DV Cost [km/s]')
ax.LineWidth = 1.08;
legend('', 'Nominal')
ylim([0, 5])

% Cobweb plot test [Not yet developed]



