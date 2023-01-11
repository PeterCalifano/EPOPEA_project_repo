%% Set Latex font
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
% This script has the aim of computing equilibrium temperature of RTG
%%
clearvars; clc; close all

% Data
% BOL Thermal power dissipated by RTG on orbiter and lander [W]
Wt_orb = 3040;
Wt_land = 1860;
sigma_SB = 5.67e-8;
% Area of the box to radiate
% Fins ?
L1_orb = 1090e-3;
L2_orb = 500e-3;
L3_orb = 500e-3;
A_orb_max = 2*L1_orb*L2_orb + L1_orb*L3_orb + 2*L3_orb*L2_orb;

L1_land = 870e-3;
L2_land = 510e-3;
L3_land = 500e-3;
A_land_max = 2*L1_land*L3_land + L1_land*L2_land + 2*L3_land*L2_land;
A_land_min = 1181238e-6;
A_orb_min = 1181238e-6;

A_orb_vec = linspace(A_orb_min,A_orb_max,100);
A_land_vec = linspace(A_land_min,A_land_max,100);
% variable: emissivity epsilon
eps_vec = linspace(0.5,1,100);

% compute equilibrium temperature 
for j = 1:100
    A_land = A_land_vec(j);
    A_orb = A_orb_vec(j);
    for k = 1:length(eps_vec)
        eps = eps_vec(k);
        T_land(j,k) = (Wt_land/(A_land*eps*sigma_SB))^(1/4);
        T_orb(j,k) = (Wt_orb/(A_orb*eps*sigma_SB))^(1/4);
    end
end

figure
surf(eps_vec, A_orb_vec, T_orb-273.15);
grid on
xlabel('$\epsilon$')
ylabel('$A\ m^2$')
zlabel('$T_{RTG}\ [°C]$')
title('Orbiter: Equilibrium temperature of RTG')

figure
surf(eps_vec, A_land_vec, T_land-273.15);
grid on
xlabel('$\epsilon$')
ylabel('$A\ m^2$')
zlabel('$T_{RTG}\ [°C]$')
title('Lander: Equilibrium temperature of RTG')

%% Heaters/ RHU required by external payload
% NAC
A_NAC = 612863e-6;
eps_NAC = 0.034; % aluminum (ASMAD)
T_NAC= -70 + 15 + 273.15;
Q_NAC = A_NAC*sigma_SB*T_NAC^4*eps_NAC;

% WAC
A_WAC = 983474e-6;
eps_WAC = 0.034; % aluminum (ASMAD)
T_WAC= - 30 + 15 + 273.15;
Q_WAC = A_WAC*sigma_SB*T_WAC^4*eps_WAC;

% TES
eps_TES = 0.034; % aluminum (ASMAD)
l1 = 130e-3; l2 = 180e-3; l3 = 180e-3;
A_TES = 23400e-6;
T_TES = 10 + 15 + 273.15;
Q_e_TES = A_TES*sigma_SB*T_TES^4*eps_TES; % heat needed to warm up TES

% Laser altimeter
l1 = 600e-3;
l2 =  400e-3; l3 = 250e-3;
A_LA = 100000e-6;
eps_LA = 0.034; % aluminum (ASMAD)
T_LA = - 20 + 15 + 273.15;
Q_LA = A_LA*sigma_SB*T_LA^4*eps_LA;

% Radar sounder
A_RS = 20341*2*(10^(-6))+21834e-6;
eps_RS = 0.034; % aluminum (ASMAD)
T_RS = - 40 + 15 + 273.15;
Q_RS = A_RS*sigma_SB*T_RS^4*eps_RS;

% compute radiators area for p/l

%% compute conductive coupling with the internal node
% COLD
% temperature of external p/l: face 3
% NAC
L2 = 4.51;
L3 = 2.7; A3 = L3*L2;
A_NAC_hot = pi*(150e-3/2)^2 ;
Q_3int3_ext_NAC = 43.4060*A_NAC_hot/A3;
T3_cold = 281.6041;
C.C_3int3ext =  2.479237200000000e+03;
T3_ext_cold = T3_cold - Q_3int3_ext_NAC/ C.C_3int3ext;
T_NAC_try = (Q_3int3_ext_NAC/eps_NAC/(A_NAC*sigma_SB))^(1/4);
Q_NAC = Q_NAC - Q_3int3_ext_NAC;

% % WAC
L2 = 4.51;
L3 = 2.7; A3 = L3*L2;
A_WAC_hot = pi*(600e-3/2)^2 ;
Q_3int3_ext_WAC = 43.4060*A_WAC_hot/A3;
T3_cold = 281.6041;
C.C_3int3ext =  2.479237200000000e+03;
T3_ext_cold = T3_cold - Q_3int3_ext_WAC/ C.C_3int3ext;
T_WAC_try = (Q_3int3_ext_WAC/eps_WAC/(A_WAC*sigma_SB))^(1/4);
Q_WAC = Q_WAC - Q_3int3_ext_WAC;

% TES
L2 = 4.51;
L3 = 2.7; A3 = L3*L2;
A_TES_hot = 23400e-6 ;
Q_3int3_ext_TES = 43.4060*A_TES_hot/A3;
T3_cold = 281.6041;
C.C_3int3ext =  2.479237200000000e+03;
T3_ext_cold = T3_cold - Q_3int3_ext_TES/ C.C_3int3ext;
T_TES_try = (Q_3int3_ext_TES/eps_TES/(A_TES*sigma_SB))^(1/4);
Q_TES = Q_e_TES - Q_3int3_ext_TES;

% laser altimeter
Q_3int3_ext_LA = 43.4060*A_LA/A3;
T3_cold = 281.6041;
C.C_3int3ext =  2.479237200000000e+03;
T3_ext_cold = T3_cold - Q_3int3_ext_LA/ C.C_3int3ext;
T_LA_try = (Q_3int3_ext_LA/eps_LA/(A_LA*sigma_SB))^(1/4);
Q_LA = Q_LA - Q_3int3_ext_LA;

Q_pl_tot = Q_NAC + Q_WAC + Q_TES + Q_RS + Q_LA;
%%
T3_ext_cold - 273.15
% hot
Q_3int3_ext_hot = -232.3214;
T3_hot =  292.3265 ;
T3_ext_hot = T3_hot - Q_3int3_ext_hot/ C.C_3int3ext;
T3_ext_hot - 273.15