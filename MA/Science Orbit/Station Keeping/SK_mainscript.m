%%%%%%%%%%%%%%%%%%%%% STATION KEEPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Load SPICE Kernels
cspice_kclear();
try
    cspice_furnsh('spice_kernels/pck00010.tpc')
    cspice_furnsh('spice_kernels/naif0012.tls')
    cspice_furnsh('spice_kernels/gm_de431.tpc')
    cspice_furnsh('spice_kernels/de440s.bsp')
    cspice_furnsh('spice_kernels/sat441.bsp')
catch
    cspice_furnsh('..\..\spice_kernels/pck00010.tpc')
    cspice_furnsh('..\..\spice_kernels/naif0012.tls')
    cspice_furnsh('..\..\spice_kernels/gm_de431.tpc')
    cspice_furnsh('..\..\spice_kernels/de440s.bsp')
end

%% DATA

% Saturn and Enceladus Data
R_Saturn = astroConstants(26);
mu_Saturn = astroConstants(16);
R_Enceladus = mean(cspice_bodvrd('602','RADII',3));
mu_Enceladus = cspice_bodvrd('602','GM',1);
J2_Saturn = 1.629061510215236e-2; % [Km^3/s^2]
J2_Enceladus = 5435.2e-6; % [Km^3/s^2]









