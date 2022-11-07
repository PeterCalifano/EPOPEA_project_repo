%% SATURN INSERTION PHASE ESTIMATIONS
clear; close all; clc

%% Load SPICE Kernels
cspice_kclear();
cspice_furnsh('..\..\spice_kernels/pck00010.tpc')
cspice_furnsh('..\..\spice_kernels/naif0012.tls')
cspice_furnsh('..\..\spice_kernels/gm_de431.tpc')
cspice_furnsh('..\..\spice_kernels/de440s.bsp')
%%
TU = 365*24*3600;
DU = 1.495978707e+8; %[Km]
date_1 = '2030-01-01 00:00:00.00 UTC';
t1 = cspice_str2et(date_1)/TU;

%% DeltaV capture estimation
% Vinf_entry = 6;
RingF_distance = 140180/R_Saturn; % [Saturn RADII]
RingG_distance(1)= 166000/R_Saturn; % [Saturn RADII]
RingG_distance(2) = 175000/R_Saturn;
id = 1;

n_points_Rp = 500;
n_points_Ra = 500;
DV_scale = 0.2:0.1:4.8;
DV_colormap = jet(length(DV_scale));

for Vinf_entry = [3, 6]
    figure(id)
    Ra_range = linspace(20, 200, n_points_Ra);
    Rp_range = linspace(RingF_distance(1), 21, n_points_Rp);
    [Ra_cap, Rp_cap] = meshgrid(Ra_range, Rp_range);

    dV_capture = nan(length(Rp_range), length(Ra_range));

    for i = 1:length(Ra_range)
        for j = 1:length(Rp_range)
            % Plot only if between F-G rings or above G ring
            %             if (Rp_range(j) > RingF_distance && Rp_range(j) < RingG_distance(1)) || Rp_range(j) > RingG_distance(2)
            dV_capture(j, i) = EstimateDVtoCapture(Vinf_entry, cspice_bodvrd('SATURN', 'GM', 1), R_Saturn*Ra_range(i), R_Saturn*Rp_range(j));
            %             end
        end
    end

    %     contour(Ra_cap, Rp_cap, dV_capture, 100, 'LineWidth', 1.05);

    colormap jet
    imagesc(Ra_range, Rp_range, dV_capture);
    hold on;
    patch([Ra_range(1), Ra_range(1), Ra_range(end), Ra_range(end)],...
        [RingG_distance(1), RingG_distance(2), RingG_distance(2), RingG_distance(1)] , [1 1 1])
    set(gca,'YDir','normal')
    view([0 0 1])
    DVcolorbar = colorbar;
    caxis([DV_scale(1), DV_scale(end)]);
    DVcolorbar.Label.String = 'DV [km/s]';
    DVcolorbar.Ticks = linspace(DV_scale(1), DV_scale(end), 15);
    axis tight
    xlabel('Apoapsis Radius [Saturn RADII]');
    ylabel('Periapsis Radius [Saturn RADII]');
%     title(['DV capture for changing capture orbit - Vinf = ', num2str(Vinf_entry), ' km/s'])
    id = id + 1;
end