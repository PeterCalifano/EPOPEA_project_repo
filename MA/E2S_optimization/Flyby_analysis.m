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

% Main backup
% load('E_EEJ_S_bestvalue1.14backup2035.mat')

% Secondary backup
% load('E_EEJ_S_bestvalue1.02backup2043.mat')

clearvars -except min_at_iter NLPoptset_local planet_seq min_pos
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
RelPos = cell(1, howmanyfb);
Eclipse_index = cell(1, howmanyfb);
time = cell(1, howmanyfb);
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

    [fb_ToF(idfb), timegrids{idfb}, xstate_cell{idfb}, Sb_cell{idfb}, SAA_cell{idfb}, RelPos{idfb}, Eclipse_index{idfb}] = EvalFlyBy(vinf_minunsplus(idfb, :), R_SOI(idfb), rp, GM_planets(idfb), Xplanets(:, idfb), bodynm);

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


    % Define directions for plot
    SunDir = -Xplanets(1:3, idfb)./norm(Xplanets(1:3, idfb));
    R_traj = xstate_cell{idfb}(:, 1:3);
    Vpdir = Xplanets(4:6, idfb)./norm(Xplanets(4:6, idfb));
    V_HelioPlus = Xplanets(4:6, idfb) + vinf_minunsplus(idfb, 4:6);
    V_HelioPlus = V_HelioPlus./norm(V_HelioPlus);

    V_HelioMinus = Xplanets(4:6, idfb) + vinf_minunsplus(idfb, 1:3);
    V_HelioMinus = V_HelioMinus./norm(V_HelioMinus);

    % Plot Sun direction and planet velocity
    quiver3(0, 0, 0, Vpdir(1), Vpdir(2), Vpdir(3), 0.3*max(vecnorm(R_traj, 2, 2)), 'Color', '#8ae222', 'LineWidth', 1.05, 'MaxHeadSize', 2);
    quiver3(0, 0, 0, SunDir(1), SunDir(2),SunDir(3), 0.3*max(vecnorm(R_traj, 2, 2)), 'Color', '#e27722', 'LineWidth', 1.05, 'MaxHeadSize', 2);

    % Plot trajectory
    idmid = (length(R_traj(:, 1))-1)/2;

    % Inbound and Entry Helio Velocity
    plot3(R_traj(1:idmid, 1), R_traj(1:idmid, 2), R_traj(1:idmid, 3), 'Color', [0.01, 0.15, 0.5], 'LineWidth', 1.05)
    quiver3(R_traj(1, 1), R_traj(1, 2), R_traj(1, 3), V_HelioMinus(1), V_HelioMinus(2), V_HelioMinus(3), 0.3*max(vecnorm(R_traj, 2, 2)), 'Color', [0.01, 0.15, 0.5], 'LineWidth', 1.05, 'MaxHeadSize', 2);
    
    % Outbound and Exit Helio Velocity
    plot3(R_traj(idmid+1:end, 1), R_traj(idmid+1:end, 2), R_traj(idmid+1:end, 3), 'Color', [0.6, 0.05, 0.01], 'LineWidth', 1.05)
    quiver3(R_traj(end, 1), R_traj(end, 2), R_traj(end, 3), V_HelioPlus(1), V_HelioPlus(2), V_HelioPlus(3), 0.3*max(vecnorm(R_traj, 2, 2)), 'Color', [0.6, 0.05, 0.01], 'LineWidth', 1.05, 'MaxHeadSize', 2);
    
    % Part of trajectory in Eclipse
    plot3(R_traj(Eclipse_index{idfb}, 1), R_traj(Eclipse_index{idfb}, 2), R_traj(Eclipse_index{idfb}, 3),'LineWidth', 1.05)

    grid minor;
    axis auto;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    xlabel('$X_{planet}$ [km]')
    ylabel('$Y_{planet}$ [km]')
    zlabel('$Z_{planet}$ [km]')

    ax.LineWidth = 1.08;
    legend(bodynm, '$V_{planet}$', 'Sun direction', 'Inbound leg', '$V_{SC}^-$', 'Outbound leg', '$V_{SC}^+$','Eclipse');
    title("Flyby " + num2str(idfb) + " at " + bodynm )
    
    % Compute Eclipse 

    % Plot angles
    i_zen = RelPos{idfb}(:,1);
    i_tan = RelPos{idfb}(:,2);
    i_tran = RelPos{idfb}(:,3);
    i_out = RelPos{idfb}(:,4);
    
    time{idfb} = [timegrids{idfb}(1:end); timegrids{idfb}(2:end) + timegrids{idfb}(end)];
    time_plot = time{idfb}/3600; % hours
    figure;
    subplot(2,2,1)
    plot(time_plot, i_zen, 'b', 'Linewidth', 1.2)
    hold on
    plot(time_plot(Eclipse_index{idfb}), i_zen(Eclipse_index{idfb}),  'Linewidth', 1.2)
    legend('$\theta_{Sun/Zenith}\ [deg]$','Eclipse')
    ylabel('$\theta_{Sun/Zenith}\ [deg]$')
    xlabel('$t\ [hours]$')
    title("Flyby " + num2str(idfb) + " at " + bodynm + "Angle between Sun and Zenith")
    grid on; grid minor

    subplot(2,2,2)
    plot(time_plot, i_out, 'g', 'Linewidth', 1.2)
    hold on
    plot(time_plot(Eclipse_index{idfb}), i_out(Eclipse_index{idfb}),  'Linewidth', 1.2)
    legend('$\theta_{Sun/Out-of-plane}\ [deg]$','Eclipse')
    ylabel('$\theta_{Sun/Out-of-plane}\ [deg]$')
    xlabel('$t\ [hours]$')
    title("Flyby " + num2str(idfb) + " at " + bodynm + "Angle out of plane")
    grid on; grid minor

    subplot(2,2,3)
    plot(time_plot, i_tan, 'g', 'Linewidth', 1.2)
     hold on
    plot(time_plot(Eclipse_index{idfb}), i_tan(Eclipse_index{idfb}),  'Linewidth', 1.2)
    legend('$\theta_{Sun/Tan}\ [deg]$','Eclipse')
    ylabel('$\theta_{Sun/Tan}\ [deg]$')
    xlabel('$t\ [hours]$')
    title("Flyby " + num2str(idfb) + " at " + bodynm + "Angle between Sun and Tangential direction")
    grid on; grid minor

    subplot(2,2,4)
    plot(time_plot, i_tran, 'g', 'Linewidth', 1.2)
     hold on
    plot(time_plot(Eclipse_index{idfb}), i_tran(Eclipse_index{idfb}),  'Linewidth', 1.2)
    legend('$\theta_{Sun/Transv}\ [deg]$','Eclipse')
    ylabel('$\theta_{Sun/Transv}\ [deg]$')
    xlabel('$t\ [hours]$')
    title("Flyby " + num2str(idfb) + " at " + bodynm + "Angle between Sun and Transversal direction")
    grid on; grid minor

end



