close all;
clear;
clc;

% Masses and volume: NSOSL, SLSL
Mass_arch = [3891; 5072];
V_estimated = 0.04*Mass_arch;

J_NSO = diag([6023.8, 11143.29, 12758.46]); % [kg m^2]
J_SO = diag([8748.36, 17399.17, 21759.28]); % [kg m^2]
J = {J_NSO; J_SO};

cm_NSO = [-0.2; 1.58; 0.01]; %[m]
cm_SO = [0.05; 1.74; 0.17]; %[m]
cm = [cm_NSO; cm_SO];

%% SRP

% Area
A_NS = [18.1, 16.705, 20.5];
A_S = [22.9, 21.875, 18];
A = [A_NS; A_S];

% Pressure center
disp_pc_cm = 3.7567;                    % [m]

% Solar Radiation Pressure
AU = [9.6 1];                         % Saturn and Earth
c = 299792458;                          % [m/s]
I = 0;                                  % incidence angle
ref = 0.5;

T_SRPmax = zeros(2,2);
for j = 1:2               % rows location: Saturn - Earth
    for id = 1:2          % cols architecture: NSO - SO
        Area = max(A(id,:));
        SolarP = 1361*(1/AU(j)^2)/c;               % Pressure
        Fsrp = SolarP*Area*(1+ref)*cos(I);
        T_SRPmax(j,id) = Fsrp*disp_pc_cm;
    end
end

fprintf('\nSRP\n')
SRP = array2table([T_SRPmax(1,:); T_SRPmax(2,:) ],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Saturn','Earth'}); 
% Display table
disp(SRP) 

%% Magnetic Disturbance

% Magnetic Dipole
class = 1e-3;           % [Am^2/kg]
D = Mass_arch.*class;   % [Am^2]

% Magnetic field of Saturn at Enceladus distance (double check and report
% reference)
B = 325*1e-9;           % [T]

T_Mmax = (D.*B)';

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
T_GGmax = [sum(T_GGMax_div(:,1)) sum(T_GGMax_div(:,2))];   % NSO - SO

fprintf('\nGRAVITY GRADIENT\n')
GG = array2table([T_GGMax_div(1,:); T_GGMax_div(2,:); T_GGmax ],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Enceladus', 'Saturn', 'Total'});
disp(GG);

%% TOTAL

% Interplanetary: SRP
Ttot_int = T_SRPmax(2,:);


% At Enceladus: GG + M +SRP
Ttot_Enc = T_SRPmax(1,:) + T_GGmax + T_Mmax;

fprintf('\nTOTAL TORQUE\n')
total = array2table([Ttot_int; Ttot_Enc],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Interplanetary', '@ Enceladus'});
disp(total);

%% Acturator Sizing

%% Reaction Wheel

% Disturbance Rejection
rotation = deg2rad(30);            %[rad]
time = 600;                        %[s]

% Slew
for i = 1:2          % NS - S
    Ji = max(diag(J{i,:}));
    T_RW_slw(i) = 4*rotation*Ji/time^2;
end

% Momentum storage
orb_period = 12*3600;

for i = 1:2
    H(i) = T_GGmax(i)*orb_period/4;
end
%% Thrusters

% Sizing for external disturbances
