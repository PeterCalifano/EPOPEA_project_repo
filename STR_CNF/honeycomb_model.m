% Model Al honeycomb structure as a uniform orthotropic material
% baseline: Al5056 
% The two face skins are activated in shear, the honeycomb core in
% compression

clear; close all; clc;

% Two face skins of full Al5056, 0.5mm each + 19mm of honeycomb core
th_face = 0.5e-3; % [m] thickness
th_core = 19e-3;

rho_Al5056 = 2640; % [kg/m^3]
rho_face = rho_Al5056;
rho_core = 130; % [kg/m^3] equivalent honeycomb density

rho_eq = (rho_core*th_core + rho_face*2*th_face)/(th_core + 2*th_face); % 20% margin considered

% Define geometry of single honeycomb cell
theta = pi/6;
d = 3e-3; % [mm] cell size
L = d/(2*cos(theta)); % [mm] single exagon side length
t = 0.01*L; % thickness of cell wall

% Compute mechanical properties
Ec = 71e9; % elastic modulus
Gc = 25.9e9; % shear modulus
nu_c = 0.33; % Poisson's ratio

Exx = Ec*(t/L)^3*cos(theta)/((1+sin(theta))*sin(theta)^2);
Eyy = Ec*(t/L)^3*(1+sin(theta))/cos(theta)^3;
nu_xy = cos(theta)^2/((1+sin(theta))*sin(theta));

Gxy = Ec*(t/L)^3*(1+sin(theta))/(2*cos(theta));

Ezz = rho_core/rho_Al5056*Ec;

nu_xz = nu_c*Exx/Ezz;
nu_yz = nu_c*Eyy/Ezz;

Gxz = Gc*(t/L)*cos(theta)/(1+sin(theta));
Gyz = Gxz;
