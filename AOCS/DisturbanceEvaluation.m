close all;
clear;
clc;

% Max lever arm length for Aerodynamic and Solar Pressure torques
SLS_Dmax = 5; % [m]
% Masses and volume: NSOSL, SLSL
Mass_arch = [3891; 5072];
V_estimated = 0.04*Mass_arch;

SLS_Vmax = 229; % [m^3]
SLS_Lmax = fzero(@(length) SLS_Vmax - pi * (SLS_Dmax/2)^2 * length, 8); % [m]

% Solar Radiation Pressure
c = 299792458; % [m/s]
SolarP = @(AU) 1361*(1/AU^2)/c;

reflectivity = 0.5;

T_SRP_max = @(AU, reflectivity, surface, length) SolarP(AU)*surface*(1+reflectivity)*(length);

distance = [9.6, 0.7]; % Saturn and Venus

for j = 1:2
    for id = 1:2

        Lmax(id) = fzero(@(length) V_estimated(id) - pi * (SLS_Dmax/2)^2 * length, 8); % [m]
        MaxLevArm(id) = Lmax(id)/2; % [m]
        Surface(id) = Lmax(id) .* SLS_Dmax;
        T_SRPmax(j, id) = T_SRP_max(distance(j), reflectivity, Surface(id), MaxLevArm(id));

    end
end

