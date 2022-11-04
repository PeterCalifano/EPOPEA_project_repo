%% EARTH-SATURN OPTIMIZATION
clear; close all; clc

%% Load SPICE Kernels
cspice_kclear();
cspice_furnsh('spice_kernels/pck00010.tpc')
cspice_furnsh('spice_kernels/naif0012.tls')
cspice_furnsh('spice_kernels/gm_de431.tpc')
cspice_furnsh('spice_kernels/de440s.bsp')
%%
list = {'Earth','Venus','Earth','Jupiter Barycenter','Saturn Barycenter'};


