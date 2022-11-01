clear;close all;clc;
% TISSERAND PLANE FOR SATURN SYSTEM
R_S = 58232;
r_a_vec = R_S*linspace(4,200);
r_p_vec = R_S*linspace(2,20);
[RA,RP] = meshgrid(r_a_vec,r_p_vec);

% Compute Tisserand parameter
C_tiss = 2./(RA+RP) + 2 * sqrt((RA+RP)/2 .* (1 - ((RA-RP)./(RA+RP)).^2));

% Conversion to Km/s
C_tiss = C_tiss/1000;

% Compute v_inf 
v_inf = sqrt(3-C_tiss);

contour(RA/R_S,RP/R_S,v_inf,10)
grid on
colorbar 
xlabel('Apoapsis [R_S]')
ylabel('Periapsis [R_S]')