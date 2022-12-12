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
l_str = 0.02;

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
epsilon_int = 0.23; %%%%%%%%%%????????????????' aluminum ????????
epsilon_ext = 0.62; % aluminized kapton see book

% MLI 
epsilon_MLI = 0.03;
alpha_MLI = 0.0029;

% Radiators
eps_rad = 0.8;
A_rad = 1*0.2;
n_rad = 4;

% RHU
n_RHU = 50;
P_RHU = 1; % 1 W is the thermal power generated by each RHU 

%%%%%%%%%%%%%%%%%%%% DEPOSIT OF DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S orbiter
%P_cold_Sorb_diss = 700;
% non sampling orbiter
%P_cold_NSorb_diss = 300; % SK
%P_cold_NSorb_plus = 374-P_cold_NSorb_diss;
% sampling orbiter
%P_cold_Sorb_plus = 866-P_cold_Sorb_diss;
%P_RTG_NSorb = 3176;
%P_RTG_Sorb = 5771;
%%%%%%%% 1. COLD CASE LANDER: ENCELADUS SURFACE %%%%%%%%%%%%%
% Power dissipated internally: POWER budget - TTMTC

% lander
P_hot_land = 513 + 20.1; % Electrical power + P dissipated by batteries
P_cold_land = 389;

% Thermal power RTG to radiate to DS / heat S/C: 
%P_RTG_land = 5191;

% Area
L = 2;
A_T = L^2; A_B = L^2; 
A_L = 4*A_T; 

% orbiter eclipses and illumination hours

% view factor
F_BT = 1/5; F_TB = 1/5; F_LT = 1/5; F_LB = 1/5;  
F_BL = 4/5; F_TL = 4/5;

%%%%%%%%% CONDUCTIVE COUPLING (C) AND RADIATIVE COUPLING (R) %%%%%%%%%%%%%%

% Sat - Saturn
% DS - Deep Space

%%% Bottom Surface (B)
R_EB1 = epsilon_Enc*A_B*epsilon_ext; % from enceladus to external surface of bottom MLI
R_B12 = A_B*epsilon_MLI; % from external to internal MLI bottom surface
R.R_EB = (1/R_EB1+1/R_B12)^(-1); % from enceladus to bottom internal surface

C_LB = 4 * k_str*(l_str*L/(L/2) + l_str*L/(L/2));
R.R_LB = A_L*F_LB*epsilon_int^2;
R.R_TB = A_T*F_TB*epsilon_int^2;

%%% Top Surface (T)
R.R_BT = A_B*F_BT*epsilon_int^2;
C_LT = 4 * k_str*(l_str*L/(L/2) + l_str*L/(L/2));
R.R_LT = A_L*F_LT*epsilon_int^2;

R_T21 = A_T * epsilon_MLI;
R_T2DS = A_T * epsilon_ext;
R.R_TDS = (1/R_T21+1/R_T2DS)^(-1);

%%% Lateral Surfaces (L)
R.R_BL = A_B*F_BL*epsilon_int^2 ;
R.R_TL = A_T*F_TL*epsilon_int^2 ;

R_L21 = (A_L - n_rad*A_rad) * epsilon_MLI;
R_L2DS = (A_L- n_rad*A_rad) * epsilon_ext;
R_LDS_rad = n_rad*A_rad * eps_rad; 
R.R_LDS = (1/R_L21+1/R_L2DS)^(-1) + R_LDS_rad;

%%% Heat Flows
Q_Sat_T = sin(pi/4)*A_T*F_sat*sigma_SB*T_Sat^4*epsilon_sat*epsilon_MLI;
Q_Sat_L = cos(pi/4)*0.5*((A_L-A_rad)*epsilon_MLI + A_rad * eps_rad)*F_sat*sigma_SB*T_Sat^4*epsilon_sat;

%%% Select powers to consider
Q_diss =  P_cold_land ;
%Q_RTG = P_RTG_land*epsilon_MLI;
Q_RTG = 0;

%%% SOLVE THE SYSTEM 
T_guess = 273*ones(3,1);
options = optimoptions('fsolve','display','iter','MaxFunctionEvaluations',50000,'Maxiterations',50000);

% Cold case
T_land_cold = fsolve(@(T) HeatBalance_Lander(T, R, T_Enc, Q_Sat_T,Q_Sat_L , Q_diss, Q_RTG, sigma_SB,C_LB,C_LT),T_guess, options);
T_land_cold-273


% Hot case
T_land_hot = fsolve(@(T) HeatBalance_Lander(T, R, T_Enc, Q_Sat_T,Q_Sat_L , Q_diss, Q_RTG, sigma_SB,C_LB,C_LT),T_guess, options);






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2. SO COLD CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7 - Antenna
% 1 - Bottom Inside
% 2 - Top Inside
% 3 - 2 Lateral with RTGs Inside
% 4 - Lateral Antenna Inside
% 5 - Lateral Instruments Inside
% 6 - Lateral Instruments Outside

clear R

L_1 = 2;
L_2 = 4;
L_3 = 4;
L = [L_1,L_2,L_3];

A_v = [L_2*L_1,L_1*L_2,L_3*L_2 + L_3*L_1,L_3*L_2,L_3*L_1];
A_tot = sum(A_v);

F_v = zeros(5,5);

for i = 1:5
    for j = 1:5
    
        F_v(i,j) = A_v(i)/(A_tot - A_v(j));

    end
end

%%% External resistences
for i = 1:5
    R_ia = A_v(i)*epsilon_ext; 
    R_ib = A_v(i)*epsilon_MLI; 
    R(i).ext = (1/R_1a + 1/R_1b)^(-1);
end

%%% Conductive Couplings
C = zeros(5,5);

for i = 1:5
    for j = 1:5
        if j ~= i
            C(i,j) = k_str*(l_str*L/(L2/2)+ l_str*L1/(L3/2)) ;
        end
    end
end

%%% Radiative Couplings inside
R_coup = zeros(5,5);
for i = 1:5
    for j = 1:5
        if i ~= j 
            R_coup(i,j) = A_v(i) * F_v(i,j) * epsilon_int^2;
        end
    end
end



%%% Heat Flows
Q_Sat_T = sin(pi/4)*A_T*F_sat*sigma_SB*T_Sat^4*epsilon_sat*epsilon_MLI;
Q_Sat_L = cos(pi/4)*0.5*((A_L-A_rad)*epsilon_MLI + A_rad * eps_rad)*F_sat*sigma_SB*T_Sat^4*epsilon_sat;

%%% Select powers to consider
Q_diss =  P_cold_land_diss + n_RHU*P_RHU;
Q_RTG = P_RTG_land*epsilon_MLI;
Q_RTG = 0;










