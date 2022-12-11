clearvars; clc; close all
% COMPUTE FLUXES

% solar fluxes
J_SUN = 1367.5;

q_sun_earth = J_SUN;
q_sun_Enc = J_SUN/9.5^2;


% Albedo 
a_earth = 0.35;
a_enc = 0.8;
a_sat = 0.499;


% View Factors 
R_earth = 6371;
R_orbit_earth = 652.4 + R_earth;
R_enc = 252.1;
R_orbit_enc_min = 20 + R_enc;
R_orbit_enc_max = 60 + R_enc;
R_sat = 58232;
R_orbit_sat = 238000;
F_earth = (R_earth/R_orbit_earth)^2;
F_enc_max = (R_enc/R_orbit_enc_min)^2;
F_enc_min = (R_enc/R_orbit_enc_max)^2;
F_sat = (R_sat/R_orbit_sat)^2;

% Radiative parameters
sigma_SB = 5.67e-8;
epsilon_Earth = 0.85;
epsilon_Enc = 1;
epsilon_sat = 1;
T_earth = 255.25;
T_Enc = 72;
T_Sat = 97;

q_albedo_earth = q_sun_earth*a_earth*F_earth;
q_albedo_sat = q_sun_Enc*a_sat*F_sat;
q_albedo_enc = q_sun_Enc*a_enc*F_enc_max; % in "hot" case

q_IR_earth = sigma_SB*epsilon_Earth*T_earth^4*F_earth;
q_IR_sat = sigma_SB*epsilon_sat*T_Sat^4*F_sat;
q_IR_enc_max = sigma_SB*epsilon_Enc*T_Enc^4*F_enc_max;
q_IR_enc_min = sigma_SB*epsilon_Enc*T_Enc^4*F_enc_min;