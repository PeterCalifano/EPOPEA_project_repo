%% Visibility of earth from Saturn's south emisphere regions
%LITERALLY ELENA'S ILLUMINATION CODE - credit to Pilo.
% SATURN USED 

clear
clc
close all

cspice_kclear()

% CHANGE WITH YOUR SCRIPT DIR: IT MUST CONTAIN ALSO THE KERNEL FOLDER
% cd('C:\Users\26ele\Desktop\POLI 2\ASMAD\PROJECT\Matlab Codes')

% kernelpool = fullfile('..\EPOPEA_metakernel.tm'); % for Windows
% cspice_furnsh(kernelpool);

%% Define Kernel
cspice_furnsh('..\spice_kernels\naif0012.tls');
cspice_furnsh ('..\spice_kernels\de440s.bsp');
cspice_furnsh ('..\spice_kernels\sat441.bsp');
cspice_furnsh ('..\spice_kernels\pck00010_mod.tpc');
cspice_furnsh ('spice_kernels\LATs.bsp');
cspice_furnsh ('spice_kernels\LATs.tf');




%% Define initial time, final time, and grid of analysis:
rot = 24*3600;
et0_obs = cspice_str2et('Jan 01 00:00:00 UTC 2030');
etf_obs = cspice_str2et('Jan 01 00:00:00 UTC 2060');
ett_obs = et0_obs:rot:etf_obs;

max_el90 = zeros(1, length(ett_obs));
max_el80 = zeros(1, length(ett_obs));
max_el60 = zeros(1, length(ett_obs));
max_el0 = zeros(1, length(ett_obs));
for i=1:length(ett_obs)
    if i == length(ett_obs)
        continue
    end
    t1 = ett_obs(i);
    t2 = ett_obs(i+1);
    tspan = t1:60:t2;

    state90 = cspice_spkpos('EARTH', tspan, 'LAT90_TOPO', 'NONE', 'LAT90');
    state80 = cspice_spkpos('EARTH', tspan, 'LAT80_TOPO', 'NONE', 'LAT80');
    state60 = cspice_spkpos('EARTH', tspan, 'LAT60_TOPO', 'NONE', 'LAT60');
    state0 = cspice_spkpos('EARTH', tspan, 'LAT0_TOPO', 'NONE', 'LAT0');
    [rho90, az90, el90] = cspice_reclat(state90);
    [rho80, az80, el80] = cspice_reclat(state80);
    [rho60, az60, el60] = cspice_reclat(state60);
    [rho0, az0, el0] = cspice_reclat(state0);
    max_el90(i) = max(rad2deg(el90));
    max_el80(i) = max(rad2deg(el80));
    max_el60(i) = max(rad2deg(el60));
    max_el0(i) = max(rad2deg(el0));
end

%% Plots
close all; clc

figure;
time = ett_obs/86400;
shift = date2jd([2000,01,01,12,00,00])-date2jd([0000,01,01,12,00,00]);
time = time + shift;
plot(time, max_el90,'color', '#0072BD', 'LineWidth', 2);
hold on
plot(time, max_el80,'color', '#D95319', 'LineWidth', 2);
plot(time, max_el60,'color', '#EDB120', 'LineWidth', 2);
plot(time, max_el0,'color', '#77AC30', 'LineWidth', 2);
grid on
grid minor
ax = gca;
ax.XTick = time(1):2*365:time(end);
datetick('x', 'yyyy mm dd', 'keeplimits', 'keepticks');
xtickangle(45);
xlim([time(1) time(end)]);
ylim([0 90]);

ylabel('Maximum Elevation [deg]')
s = legend('0° long, -90° lat', '0° long, -80° lat', '0° long, -60° lat', '0° long, 0° lat');
s.Location = 'best';
