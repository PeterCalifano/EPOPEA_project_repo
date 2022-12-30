close all
clear
clc

DefaultFontSize = 16;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

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

% Normalization coefficients
% TU = 24*3600; % 1 day
% DU = 1.495978707e+8; % 1 AU

%% Load solution workspace
% Baseline
load('E_EEJ_S_bestvalue1.016.mat')

clearvars -except min_at_iter NLPoptset_local planet_seq min_pos
% Main backup
% load('E_EEJ_S_bestvalue1.14backup2035.mat')

% Secondary backup
% load('E_EEJ_S_bestvalue1.02backup2043.mat')

%% Compute trajectory
R_Saturn = astroConstants(26); % [km]
Ra_target = 200*R_Saturn;
Rp_target = 2.55*R_Saturn;

% Extract best solution
[~, index_iter] = min(min_at_iter(min_at_iter > 0));
index_pos = min_pos(index_iter);
Xsol = NLPoptset_local(index_pos, :, index_iter);

% Evaluate solution
[DVsum, data] = objfun_EStransfer_analysis(Xsol, planet_seq, Ra_target, Rp_target);

Xplanets = data.Xplanets;
EventEpochsJD = data.EventEpochsJD;
vinf_minunsplus = data.vinf_minusplus;
DV_breakdown = data.DVs;

%% Evaluate flybys
howmanyfb = numel(EventEpochsJD);
fbplanets = planet_seq(2:end-1);

R_SOI = zeros(1, howmanyfb);
GM_planets = zeros(1, howmanyfb);
MeanR_planets = zeros(1, howmanyfb);

GM_Sun = cspice_bodvrd('SUN', 'GM', 1);
fb_ToF = zeros(1, howmanyfb); % [hours]
timegrids = cell(1, howmanyfb);
xstate_cell = cell(1, howmanyfb);
Sb_cell = cell(1, howmanyfb);
SAA_cell = cell(1, howmanyfb);

for idfb = 1:howmanyfb
    % Determine Planet properties
    switch fbplanets(idfb)
        case 1
            bodynm = 'Mercury';
        case 2
            bodynm = 'Venus';
        case 3
            bodynm = 'Earth';
        case 4
            bodynm = 'Mars';
        case 5
            bodynm = 'Jupiter';
        case 6
            bodynm = 'Saturn';
    end

    GM_planets(idfb) = cspice_bodvrd(bodynm, 'GM', 1);
    MeanR_planets(idfb) = mean(cspice_bodvrd(bodynm, 'RADII', 3));
    R_SOI(idfb) = norm(Xplanets(1:3, idfb)) * (GM_planets(idfb)/GM_Sun)^(2/5);

    rp = Xsol(6 + 2*howmanyfb + idfb).*MeanR_planets(idfb); 

    [fb_ToF(idfb), timegrids{idfb}, xstate_cell{idfb}, Sb_cell{idfb}, SAA_cell{idfb}] = EvalFlyBy(vinf_minunsplus(idfb, :), R_SOI(idfb), rp, GM_planets(idfb), Xplanets(:, idfb));

    figure('WindowState', 'maximized');
    hold on;
    try
        opts.Units = 'km';
        opts.RefPlane = 'equatorial';
        planet = planet3D(bodynm, opts);
        hold on;
    catch
        plot3(0, 0, 0, 'ko','MarkerFaceColor', 'k', 'DisplayName', bodynm);
        hold on;
    end

    SunDir = Xplanets(1:3, idfb)./norm(Xplanets(1:3, idfb));

    R_traj = xstate_cell{idfb}(:, 1:3);
    plot3(R_traj(:, 1), R_traj(:, 2), R_traj(:, 3), 'Color', [0.01, 0.15, 0.5], 'LineWidth', 1.05)

    quiver3(0, 0, 0, SunDir(1), SunDir(2), SunDir(3), 0.3*max(vecnorm(R_traj, 2, 2)), 'Color', '#e27722', 'LineWidth', 1.05, 'MaxHeadSize', 2);

    grid minor;
    axis auto;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    xlabel('$X_{planet}$ [km]')
    ylabel('$Y_{planet}$ [km]')
    zlabel('$Z_{planet}$ [km]')

    ax.LineWidth = 1.08;
    legend(bodynm, 'Trajectory', 'Sun direction');

end







