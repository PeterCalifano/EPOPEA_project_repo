clc

% kernelpool = '..\EPOPEA_project_repo\MA\EPOPEA_metakernel.tm';
kernelpool = fullfile('..\EPOPEA_metakernel.tm'); % for Windows
% kernelpool = fullfile('..\EPOPEA_metakernelUNIX.tm'); % for UNIX/MAC

% kernelpool = '..\EPOPEA_metakernel.tm'
cspice_kclear();

cspice_furnsh(kernelpool);

et0 = cspice_str2et('2022-03-10 11:00 UTC');
cspice_spkezr('EARTH', et0, 'ECLIPJ2000', 'NONE', 'SSB')






%% B-plane test

% v_inf_minus = [15, 2, 0];

% Build B-plane for given input target 
% Direction of INCOMING asymptote in Planet RF
S_hat = v_inf_minus./norm(v_inf_minus); 
N_hat = [0, 0, 1]; % Direction of North Pole of the Planet (can it be taken as direction of its angular momentum?)

% Direction in the Equatorial Plane of the Planet
T_hat = cross(S_hat, N_hat)./norm(cross(S_hat, N_hat));

R_hat = cross(S_hat, T_hat)./norm(cross(S_hat, T_hat));

% B = [2, 4, 0.6];
Bt = dot(B, T_hat);
Br = dot(B, R_hat);

planetname = 'saturn';

% Plot B-plane and (Bt, Br) components
figure;
hold on;
plot(0, 0, 'Marker', '.', 'MarkerSize', 15, 'Color', PlanetColour(planetname));

plot(Bt, Br, 'Marker', 'x', 'LineWidth', 1.08, 'Color', '#dd4411');


% Plot options
axis ij % Reverse y axis direction (increase from top to bottom)
yline(0,'LineWidth', 1.02)
xline(0,'LineWidth', 1.02)
hold off;

ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
grid minor

xlabel('T [km]');
ylabel('R [km]')
axis auto




