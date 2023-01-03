%% Set Latex font
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
%%
clearvars; clc; close all

% Solar Flux
AU_earth = 1; 
AU_enc = 9.5;
q_sun_enc = 1367.5/AU_enc^2;
q_sun_earth = 1367.5;

% Albedo 
a_earth = 0.35;
a_enc = 0.8;
a_sat = 0.499;

% View Factors 
R_earth = 6371;
R_orbit_earth = 670 + R_earth;
R_enc = 252.1;
R_orbit_enc_min = 20 + R_enc;
R_orbit_enc_max = 1.05e+03 + R_enc;
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

%%% Spacecraft %%%

%%% Structure properties (Al-5056-O)
k_str = 117;
l_str = 0.02;

% Structure
epsilon_int = 0.23; %%%%%%%%%%????????????????' aluminum ???????? -> to change, sentitivity analysis

% MLI 
epsilon_MLI = 0.03;
% alpha_MLI = 0.14; % max. typ = 0.12 
alpha_MLI = 0.12;

% Antenna painted white
epsilon_ant = 0.90;
alpha_ant = 0.18; % do sensitivity analysis

% Thermal straps
k_TS = 398; 
L_TS = 120e-3; %?????

% Radiators
% changed from 0.8: if louvers open ~0.7
eps_louv_closed = 0.14;
eps_louv_open = 0.7;
eps_rad = eps_louv_open; % if open louvers

% eps_rad = 0.8;  % if ideal radiators and not louvers
alpha_louv_closed = 0.062; 
alpha_louv_open = 0.269; % (worst case EOL)
A_rad_one = 1 * 0.2; % area one radiator
n_rad = 4; %% can be changed
A_rad_tot = A_rad_one*n_rad;% total area of radiators
k_rad = 1;

% Areas
L = 1.5; L1 = L; L2 = L; L3 = L;
A1 = L^2;
A2 = L^2;
A3 = L^2;
A4 = L^2;
A5_tot = L^2;
A6 = L^2;

A5 = A5_tot-A_rad_tot;

% RHU or heaters
n_RHU = 0; % can be changed ! 
P_RHU = 1; % 1 W is the thermal power generated by each RHU 

% Power
P_budget_hot = 0;
P_input_TMTC_h = 0;
P_diss_TMTC_h = 0;
% add batteries ...
We =  456;
Wt = 3040;
eff_shunt = 0.6;
Q_hot = We*eff_shunt;
Q_hot  = 0;
P_shunt = Q_hot - We; % check power can be dissipated though a shunt
% view factor
% Surface 1
F12 = VF_PerpRec(L3,L2,L1);
F21 = F12*A1/A2;
F13 = VF_PerpRec(L3,L1,L3);
F31 = F13*A1/A3;
F14 = VF_PerpRec(L3,L2,L1);
F41 = VF_PerpRec(L2,L3,L1);
F15_tot = F13;
F51_tot = F15_tot*A5_tot/A1;
F1rad = F15_tot*A_rad_tot/A5_tot;
Frad1 = F1rad*A1/A_rad_tot; 
F15 = F15_tot*A5/A5_tot;
F51 = F15*A1/A5;
F16 = VF_ParallelEqualRec(L3,L1,L2);
F61 = F16;

% Surface 2 
F23 = VF_PerpRec(L2,L3,L1);
F32 = F23*A2/A3;
F24 = VF_ParallelEqualRec(L2,L1,L3);
F42 = F24*A2/A4;
F25_tot = VF_PerpRec(L2,L3,L1);
F52_tot = F25_tot*A5_tot/A2;
F2rad = F25_tot*A_rad_tot/A5_tot;
Frad2 = F2rad*A2/A_rad_tot; 
F25 = F25_tot*A5/A5_tot;
F52 = F25*A2/A5;
F26 = VF_PerpRec(L2,L3,L1);
F62 = F26*A2/A6;

% Surface 3
F34 = F32;
F43 = F23;
F35_tot = VF_ParallelEqualRec(L1,L2,L3);
F53_tot = F35_tot*A5_tot/A3;
F3rad = F35_tot*A_rad_tot/A5_tot;
Frad3 = F3rad*A3/A_rad_tot; 
F35 = F35_tot*A5/A5_tot;
F53 = F35*A3/A5;
F36 = F31;
F63 = F13;

% Surface 4
F45_tot = F43;
F54_tot = F45_tot*A5_tot/A4;
F4rad = F45_tot*A_rad_tot/A5_tot;
Frad4 = F4rad*A4/A_rad_tot; 
F45 = F45_tot*A5/A5_tot;
F54 = F45*A4/A5;
F46 = F41; 
F64 = F14;

% Surface 6
F65_tot = F63;
F56_tot = F65_tot*A5_tot/A6;
F6rad = F65_tot*A_rad_tot/A5_tot;
Frad6 = F6rad*A6/A_rad_tot; 
F65 = F65_tot*A5/A5_tot;
F56 = F65*A6/A5;

% radiative coupling
R.R_12 = sigma_SB * A1 * epsilon_int^2 * F12;
R.R_1rad = sigma_SB * A1 * epsilon_int * eps_louv_open * F1rad;
R.R_13 = sigma_SB * A1 * epsilon_int^2 * F13;
R.R_14 = sigma_SB * A1 * epsilon_int^2 * F14;
R.R_15 = sigma_SB * A1 * epsilon_int^2 * F15;
R.R_16 = sigma_SB * A1 * epsilon_int^2 * F16;
R.R_23 = sigma_SB * A2 * epsilon_int^2 * F23;
R.R_24 = sigma_SB * A2 * epsilon_int^2 * F24;
R.R_25 = sigma_SB * A2 * epsilon_int^2 * F25;
R.R_26 = sigma_SB * A2 * epsilon_int^2 * F26;
R.R_rad3 = sigma_SB * A_rad_tot * epsilon_int * eps_louv_open * Frad3; % ?
R.R_rad4 = sigma_SB * A_rad_tot * epsilon_int * eps_louv_open * Frad4; % ?
R.R_rad2 = sigma_SB * A_rad_tot * epsilon_int * eps_louv_open * Frad2; % ?
R.R_rad6 = sigma_SB * A_rad_tot * epsilon_int * eps_louv_open * Frad6; % ?
R.R_rad5 = 0;
R.R_34 = sigma_SB * A3 * epsilon_int^2 * F34;
R.R_35 = sigma_SB * A3 * epsilon_int^2 * F35;
R.R_36 = sigma_SB * A3 * epsilon_int^2 * F36;
R.R_45 = sigma_SB * A4 * epsilon_int^2 * F45;
R.R_46 = sigma_SB * A4 * epsilon_int^2 * F46;
R.R_56 = sigma_SB * A5 * epsilon_int^2 * F56;
R.R_10 = sigma_SB * A1 * epsilon_MLI;
R.R_20 = sigma_SB * A2 * epsilon_MLI;
R.R_30 = sigma_SB * A3 * epsilon_MLI;
R.R_40 = sigma_SB * A4 * epsilon_MLI;
R.R_50 = sigma_SB * A5 * epsilon_MLI;
R.R_60 = sigma_SB * A6 * epsilon_MLI;
R.R_rad0 = sigma_SB * A_rad_tot * eps_rad;
% Conductive Coupling
C.C_12 = k_str*l_str*L1*(1/(L2/2)+1/(L3/2));
C.C_13 = k_str*l_str*L2*(1/(L1/2)+1/(L3/2));
C.C_14 = k_str*l_str*L1*(1/(L2/2)+1/(L3/2));
C.C_15 = k_str*l_str*L2*(1/(L1/2)+1/(L3/2));
C.C_16 = 0;
C.C_23 = k_str*l_str*L3*(1/(L1/2)+1/(L2/2));
C.C_25 = k_str*l_str*L3*(1/(L1/2)+1/(L2/2));
C.C_24 = 0;
C.C_26 = C.C_12;
C.C_34 = C.C_23;
C.C_35 = 0;
C.C_36 = C.C_13;
C.C_45 = C.C_34;
C.C_46 = C.C_14;
C.C_56 = C.C_13;

C.C_5rad = k_str*(l_str*L3/(L2/2));

% to tune:
C.C_1rad = 10;
C.C_2rad = 10;
C.C_3rad = 0;    % opposite surface
C.C_4rad = 10;
C.C_6rad = 10;
C.C_5rad = C.C_5rad + 100;

% External fluxes 
q_Sun = q_sun_earth;
q_alb = q_Sun*F_earth*a_earth;
q_Earth = F_earth*sigma_SB*T_earth^4*epsilon_Earth;

% Angles with external fluxes

theta_6Sun = 24.7 * pi/180;
theta_3Sun = 90 - theta_6Sun;

Q_ext_hot = zeros(8,1);
% HOT CASE 1 : cameras towards Earth, face 5 sees the Sun
% HOT CASE 2: cameras towards Earth, antenna towards Sun
% HOT CASE 3: HGA (1) towards Earth, RTGs (6) towards Sun
hot_case = 3;

switch hot_case
    case 1
        theta_5Sun = 15 * pi/180;
        Q_ext_hot(3) = q_Earth*epsilon_MLI*A3 + q_alb*A3*alpha_MLI ;
        Q_ext_hot(5) = q_Sun * A5 * alpha_MLI* cos(theta_5Sun); 
        Q_ext_hot(6) = q_Sun * A6 * alpha_MLI* cos(theta_6Sun);
    case 2
        theta_1Sun = 0;
        Q_ext_hot(1) = q_Sun * A1 * alpha_MLI* cos(theta_1Sun);
        Q_ext_hot(3) = q_Earth*epsilon_MLI*A3 + q_alb*A3*alpha_MLI ;
        Q_ext_hot(7) = q_Sun * A_rad_tot * alpha_louv_open* cos(theta_1Sun);
    case 3
        Q_ext_hot(6) = q_Sun * A6 * alpha_MLI* cos(theta_6Sun);
        Q_ext_hot(3) = q_Sun * A3 * alpha_MLI* cos(theta_3Sun);
        Q_ext_hot(1) = q_Earth * epsilon_MLI * A1 + q_alb * A1 * alpha_MLI ;
end
% Initial condition
T0 = 293;

% solve
%%% Internal dissipation power
Q_diss_hot =  Q_hot;

%%% SOLVE THE SYSTEM 
Clamped = 1;
T_guess = 273*ones(7,1);
options = optimoptions('fsolve','display','iter','MaxFunctionEvaluations',50000,'Maxiterations',50000);
T_land_hot = fsolve(@(T) HeatBalance_Lander_new(T, R, C, Q_ext_hot , Q_diss_hot, Clamped), T_guess, options);

fprintf(['1 ',num2str(T_land_hot(1)-273),' Celsius\n'])
fprintf(['2 ',num2str(T_land_hot(2)-273),' Celsius\n'])
fprintf(['3 ',num2str(T_land_hot(3)-273),' Celsius\n'])
fprintf(['4 ',num2str(T_land_hot(4)-273),' Celsius\n'])
fprintf(['5 ',num2str(T_land_hot(5)-273),' Celsius\n'])
fprintf(['6 ',num2str(T_land_hot(6)-273),' Celsius\n'])
fprintf(['rad ',num2str(T_land_hot(7)-273),' Celsius\n'])
% add mass and specific heat for transient

%% Cold case
% Power
P_budget_cold = 277;
P_input_TMTC_cold = 20;       % ask Antoine
P_diss_TMTC_cold =18.86;      % ask Antoine
% add batteries ...
We = 389;
Wt = 2600;
Q_cold = P_budget_cold-P_input_TMTC_cold+P_diss_TMTC_cold;

Q_diss_cold = Q_cold;

% close louvers and compute again thermal couplings
% perc_open_rad = 0.5; % percentage of open radiators
% eps_rad = perc_open_rad * eps_louv_open + (1-perc_open_rad) * eps_louv_closed;
eps_rad = eps_louv_closed;
% radiative coupling
R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;

% External fluxes
% IR Heat fluxes for Saturn and Enceladus
q_Sat = F_sat*sigma_SB*T_Sat^4*epsilon_sat;
q_Enc_orbit = F_enc_min * sigma_SB *T_Enc^4 *epsilon_Enc;
q_Enc_ground = 1 * sigma_SB *T_Enc^4 *epsilon_Enc;

% HYP: nadir pointing
% face 3 with cameras points Enceladus and receive IR from Enceladus
% Saturn ? HYP: face 4 (CHANGE!)
P_added_RHU = 0;

theta_3Enc = 0;
theta_6Sat = deg2rad(20); % CHANGE!
theta_4Sat = deg2rad(70); % CHANGE!

cold_case = 1; % orbit / during landing
% cold_case = 2; % on ground
Q_ext_cold = zeros(7,1);

switch cold_case
    case 1
    theta_6Sat = 0; % CHANGE!
    Clamped = 0;
    Q_ext_cold(3) = q_Enc_orbit * A3 * epsilon_MLI * cos(theta_3Enc);
    Q_ext_cold(6) = q_Sat * A6 * epsilon_MLI*cos(theta_6Sat); 
    case 2
    theta_3Enc = deg2rad(20);
    theta_2Enc = 0;
    Clamped = 0;
    Q_ext_cold(2) = q_Enc_ground * A2 * epsilon_MLI * cos(theta_2Enc);
    Q_ext_cold(4) = q_Sat * A4*epsilon_MLI*cos(theta_4Sat);
    Q_ext_cold(6) = q_Sat * A6*epsilon_MLI*cos(theta_6Sat);
end

%%% SOLVE THE SYSTEM 
T_guess = 273*ones(7,1);
options = optimoptions('fsolve','display','iter','MaxFunctionEvaluations',50000,'Maxiterations',50000);
T_land_cold = fsolve(@(T) HeatBalance_Lander_new(T, R, C, Q_ext_cold , Q_diss_cold, Clamped), T_guess, options);

fprintf(['1 ',num2str(T_land_cold(1)-273),' Celsius\n'])
fprintf(['2 ',num2str(T_land_cold(2)-273),' Celsius\n'])
fprintf(['3 ',num2str(T_land_cold(3)-273),' Celsius\n'])
fprintf(['4 ',num2str(T_land_cold(4)-273),' Celsius\n'])
fprintf(['5 ',num2str(T_land_cold(5)-273),' Celsius\n'])
fprintf(['6 ',num2str(T_land_cold(6)-273),' Celsius\n'])
fprintf(['rad ',num2str(T_land_cold(7)-273),' Celsius\n'])

% compute the temperature of surface 2 external (where cameras are):
% C.C_2int2ext = k_str*A2/l_str;
% R.R_2int2ext = sigma_SB*A2 * epsilon_MLI;
% Q_20 = R.R_20*((T_land_cold(2))^4 - 0);
% Q_2int2_ext = Q_20;
% T2_ext = (T_land_cold(2)) - Q_2int2_ext / C.C_2int2ext; % not sure about this
% fprintf(['2 ext  ',num2str(T2_ext-273),' Celsius\n'])

% Conduction between surfaces
% To do:
% 1) add conduction 
%      - between surfaces
%      - to radiators
% 2) change internal power
% 3) surface in contact with lander?
% 4) Add internal nodes
% sensitivity analysis, check properties
% T of RTG


