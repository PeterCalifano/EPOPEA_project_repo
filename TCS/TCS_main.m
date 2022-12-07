clearvars; clc; close all
%%% DATA environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Solar Flux
AU_earth = 1; 
AU_enc = 9.5;
q_sun_enc = 1367.5/AU_enc^2;
q_sun_earth = 1367.5;

%%% Albedo 
a_earth = 0.35;
a_enc = 0.8;
a_sat = 0.499;
R_earth = 6371;
R_orbit_earth = 470 + R_earth;
R_enc = 252.1;
R_orbit_enc_min = 20 + R_enc;
R_orbit_enc_max = 20 + R_enc;
R_sat = 58232;
R_orbit_sat = 238000;

% Structure properties (Al-5056-O)
k_str = 117;
l_str = 0.02; % TBC

% Radiators
eps_rad = 1;
A_rad = 0.5*0.5;

%%% View Factors 
F_earth = (R_earth/R_orbit_earth)^2;
F_enc_max = (R_enc/R_orbit_enc_min)^2;
F_enc_min = (R_enc/R_orbit_enc_max)^2;
F_sat = (R_sat/R_orbit_sat)^2;

%%% Radiative parameters
sigma_SB = 5.67e-8;
epsilon_Earth = 0.85;
epsilon_Enc = 1;
epsilon_sat = 1;
T_earth = 255.25;
T_Enc = 72;
T_Sat = 97;

% MLI 
epsilon_MLI = 0.03;
alpha_MLI = 0.0029;

epsilon_int = 0.23; %%%%%%%%%%????????????????' aluminum ????????
epsilon_ext = 0.62; % aluminized kapton see book

%epsilon_ext = 1;


%%%%%%%% 1. COLD CASE ENCELADUS ORBITER AND LANDER SEPARATELY %%%%%%%%%%%%%
% Power dissipated internally: POWER budget - TTMTC

% non sampling orbiter
P_cold_NSorb_diss = 300; % SK
% lander
P_cold_land_diss = 731.77-201.6;
% S orbiter
P_cold_Sorb_diss = 700;

% surplus of power: we can choose it warm up S/C or to radiate
P_cold_NSorb_plus = 374-P_cold_NSorb_diss;
% lander
P_cold_land_plus = 778-P_cold_land_diss;
% sampling orbiter
P_cold_Sorb_plus = 866-P_cold_Sorb_diss;
% Thermal power RTG to radiate to DS / heat S/C: 
P_RTG_NSorb = 3176;
P_RTG_land = 5191;
P_RTG_Sorb = 5771;

% Area
L = 2;
A_T = L^2; A_B = L^2; 
A_L = 4*A_T; 

% lander always Saturn 
% lander tot hours Sun tot our not illuminated

% orbiter eclipses and illumination hours

% view factor
F_BT = 1/5; F_TB = 1/5; F_LT = 1/5; F_LB = 1/5;  
F_BL = 4/5; F_TL = 4/5;

% RADIATIVE COUPLING
% bottom surface
R_EB1 = epsilon_Enc*A_B*epsilon_ext; % from enceladus to external surface of bottom MLI
R_B12 = A_B*epsilon_MLI; % from external to internal MLI bottom surface
R.R_EB = (1/R_EB1+1/R_B12)^(-1); % from enceladus to bottom internal surface

C_LB = 4 * k_str*(l_str*L/(L/2) + l_str*L/(L/2));
R.R_LB = A_L*F_LB*epsilon_int^2;
R.R_TB = A_T*F_TB*epsilon_int^2;

%top
R.R_BT = A_B*F_BT*epsilon_int^2;
C_LT = 4 * k_str*(l_str*L/(L/2) + l_str*L/(L/2));
R.R_LT = A_L*F_LT*epsilon_int^2;

R_T21 = (A_T - A_rad) * epsilon_MLI;
R_T2DS = (A_T - A_rad) * epsilon_ext;
R_TDS_rad = A_rad * eps_rad;  
R.R_TDS = (1/R_T21+1/R_T2DS)^(-1) + R_TDS_rad;

% lateral
R.R_BL = A_B*F_BL*epsilon_int^2 ;
R.R_TL = A_T*F_TL*epsilon_int^2 ;

R_L21 = A_L*epsilon_MLI;
R_L2DS = A_L*epsilon_ext;
R.R_LDS = (1/R_L21+1/R_L2DS)^(-1);

% heat flow
Q_Sat_T = A_T*F_sat*sigma_SB*T_Sat^4*epsilon_sat*epsilon_MLI;

Q_Sat_L = A_L*F_sat*sigma_SB*T_Sat^4*epsilon_sat*epsilon_MLI;

% first guess
T_guess = 273*ones(3,1);
Q_diss =  P_cold_land_diss;
Q_RTG = P_RTG_land*epsilon_MLI;
Q_RTG = 0;
options = optimoptions('fsolve','display','iter','MaxFunctionEvaluations',50000,'Maxiterations',50000);
T_nodes = fsolve(@(T) HeatBalance(T, R, T_Enc, Q_Sat_T,Q_Sat_L , Q_diss, Q_RTG, sigma_SB,C_LB,C_LT),T_guess, options);

T_nodes-273

%%
function [balances] = HeatBalance(T, R, T_E, Q_Sat_T,Q_Sat_L , Q_diss, Q_RTG, sigma,C_LB,C_LT)

T_B = T(1);
T_L = T(2);
T_T = T(3);

Cond_LB = C_LB/(4*sigma*0.5^3*(T_L + T_B)^3);
Cond_LT = C_LT/(4*sigma*0.5^3*(T_L + T_T)^3);

% balance bottom
Q_EB = R.R_EB*sigma*(T_E^4-T_B^4);
Q_LB = (R.R_LB + Cond_LB)*sigma*(T_L^4-T_B^4);
Q_TB = R.R_TB*sigma*(T_T^4-T_B^4);

% balance lateral surfaces
Q_BL = (R.R_BL + Cond_LB) *sigma*(T_B^4-T_L^4);
Q_TL = (R.R_TL + Cond_LT) *sigma*(T_T^4-T_L^4);
Q_LDS = R.R_LDS*sigma*T_L^4;

% balance top surface
Q_BT = R.R_BT*sigma*(T_B^4-T_T^4);
Q_LT = (R.R_LT + Cond_LT)*sigma*(T_L^4-T_T^4);
Q_TDS = R.R_TDS*sigma*T_T^4;

balances(1) = Q_EB+Q_LB+Q_TB;
balances(2) = Q_BL+Q_TL-Q_LDS+Q_Sat_L+Q_diss+Q_RTG;
balances(3) = Q_BT+Q_LT-Q_TDS+Q_Sat_T; 

end
