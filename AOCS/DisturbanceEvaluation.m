close all;
clear;
clc;

% Masses and volume: NSOSL, SLSL
Mass_arch = [3891; 5072];
V_estimated = 0.04*Mass_arch;

J_NSO = diag([5177.71, 9373.88, 11328.19]); % [kg m^2]
J_SO = diag([7819.9 14791.83 19333.73]); % [kg m^2]

%% SRP

% Max lever arm length for Aerodynamic and Solar Pressure torques
SLS_Dmax = 5; % [m]

SLS_Vmax = 229; % [m^3]
SLS_Lmax = fzero(@(length) SLS_Vmax - pi * (SLS_Dmax/2)^2 * length, 8); % [m]

% Solar Radiation Pressure
c = 299792458; % [m/s]
SolarP = @(AU) 1361*(1/AU^2)/c;

reflectivity = 0.5;

T_SRP_max = @(AU, reflectivity, surface, length) SolarP(AU)*surface*(1+reflectivity)*(length);

distance = [9.6, 0.7]; % Saturn and Venus

for j = 1:2               % rows location: Saturn - Venus
    for id = 1:2          % cols architecture: NSO - SO

        Lmax(id) = fzero(@(length) V_estimated(id) - pi * (SLS_Dmax/2)^2 * length, 8); % [m]
        MaxLevArm(id) = Lmax(id)/2; % [m]
        Surface(id) = Lmax(id) .* SLS_Dmax;
        T_SRPmax(j, id) = T_SRP_max(distance(j), reflectivity, Surface(id), MaxLevArm(id));

    end
end

fprintf('\nSRP\n')
SRP = array2table([T_SRPmax(1,:); T_SRPmax(2,:) ],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Saturn','Venus'}); 
% Display table
disp(SRP) 

%% Magnetic Disturbance

% Magnetic Dipole
class = 1e-3;           % [Am^2/kg]
D = Mass_arch.*class;   % [Am^2]

% Magnetic field of Saturn at Enceladus distance (double check and report
% reference)
B = 325*1e-9;           % [T]

T_Mmax = D.*B;

fprintf('\nMAGNETIC TORQUE\n')
M = array2table([T_Mmax(1) T_Mmax(2) ],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Enceladus'}); 
disp(M)

%% Gravity Gradient

muS = 3.7931187*1e16;            %[m^3/s^2]
muE = (6.6743e-11)*(1.0802e20);  %[m^3/s^2]
mu = [muE muS];

% Distance: minimum for Enceladus, mean from Saturn
min_d_E = 29000+251500;                 %[m]
mean_d_S = 238000000;            %[m] 
dist = [min_d_E mean_d_S];

% Minimum and maximum inertia
Imin = [min(diag(J_NSO)) min(diag(J_SO))];
Imax = [max(diag(J_NSO)) max(diag(J_SO))];

% Tilting angle wrt Nadir: worst case
theta = deg2rad(45);

T_GGMax_div = zeros(2,2);
for p = 1:2         %rows: Enceladus, Saturn
    for a = 1:2     %cols: NSO, SO
        T_GGMax_div(p,a) = 3*mu(p)/(2*dist(p)^3) * abs(Imin(a)-Imax(a)) * sin(2*theta);
    end
end
T_GGMax = [sum(T_GGMax_div(:,1)) sum(T_GGMax_div(:,2))];   % NSO - SO

fprintf('\nGRAVITY GRADIENT\n')
GG = array2table([T_GGMax_div(1,:); T_GGMax_div(2,:); T_GGMax ],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Enceladus', 'Saturn', 'Total'});
disp(GG);