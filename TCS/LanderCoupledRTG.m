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
% l_str = 0.02;
l_str = 0.003;
l_str_tot = 10e-3;
k_honeycomb = 2.036; % sensitivity analysis: 0.852
% Structure
% epsilon_int = 0.23; %%%%%%%%%%????????????????' aluminum ???????? -> to change, sentitivity analysis
% epsilon_int = 0.874; % black paint to maximize exchange (you can change it)
% epsilon_int = 0.034; 
epsilon_5 =  0.874; % black paint to maximize exchange
epsilon_6 = 0.874; % black paint to maximize exchange
epsilon_1 = 0.874; 
epsilon_2 = 0.874; 
epsilon_3 = 0.874; 
epsilon_4 = 0.874;
% epsilon_8 = 0.874;
epsilon_8 = 0.034;
% epsilon_5 = 0.034;
% epsilon_6 = 0.034;
% epsilon_8 = 0.034;
% epsilon_1 = 0.034;
% epsilon_2 = 0.034;
% epsilon_3 = 0.034;
% epsilon_4 = 0.034;

% MLI 
epsilon_MLI = 0.03;
% alpha_MLI = 0.08 (Kapton, silvered, aluminum oxide coated, 1 mil)
%           = 0.11 Kapton, aluminized, silicon oxide coated, 1 mil
%           = 0.12 typical
alpha_MLI = 0.11;

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
A_rad_one = 1 *530*390e-6; % area one radiator
n_rad = 3.5; %% can be changed
A_rad_tot = A_rad_one*n_rad;% total area of radiators
k_rad = 1;

% Areas
L1 = 1.5; L3 = 1.5; L2 = 1.8;
A1 = L1*L2;
A2 = L1*L3;
A3 = L2*L3;
A4 = A2;
A5 = A3;
A6_tot = A1;
A6_int = A1;
A8 = 912319e-6; % from CAD: base of pl
A6_ext = A6_tot - A_rad_tot;


% RHU or heaters
n_RHU = 0; % can be changed ! 
P_RHU = 1; % 1 W is the thermal power generated by each RHU 

% Power
P_budget_hot = 0;
P_input_TMTC_h = 0;
P_diss_TMTC_h = 0;
% add batteries ...
We = 279;
Wt = 1860;
% eff_shunt = 0.6;
% Q_hot = We*eff_shunt;
Q_hot  = 0;
P_shunt = - Q_hot + We; % check power can be dissipated though a shunt
% view factor
% Surface 1
H8 = 830e-3; 
L = L1;
A8 = L1*L3;
F12_a = VF_PerpRec(L3,H8,L1);
F13_a = VF_PerpRec(L3,L1,H8);
F15_a = F13_a;
F16_a = VF_ParallelEqualRec(L3,L1,H8);
F18_a = VF_PerpRec(L3,H8,L1)*H8/L2;
F12 = (F12_a + F13_a + F15_a + F16_a)*H8/L2;
F13 = VF_PerpRec(L3,L1,L2-H8)*(L2-H8)/L2;
F14 = VF_PerpRec(L3,L2-H8,L1)*(L2-H8)/L2;
F15 = F13*(L2-H8)/L2;
F16 = VF_ParallelEqualRec(L3,L1,L2-H8)*(L2-H8)/L2;
F18_b = VF_PerpRec(L3,L2-H8,L1)*(L2-H8)/L2;
F18 = F18_a + F18_b;
F21 = F12_a*H8/L2;
F31 = F13;
F41 = F14;
F51 = F15;
F61 = F16;
F81 = F18*A1/A8;

% Surface 2
F23 = F21;
F32 = F12;
F24 = 0;
F42 = 0;
F25 = F21;
F52 = F12;
F26 = F21;
F62 = F12;
F28 = VF_ParallelEqualRec(H8,L3,L1);
F82 = F28;

% Surface 3
F34 = F14;
F43 = F41;
F35 = F16;
F53 = F61;
F36 = F13;
F63 = F31;
F38 = F18;
F83 = F81;

% Surface 4
F45 = F43;
F54 = F34;
F46 = F45;
F64 = F54;
F48 = VF_ParallelEqualRec(L2-H8,L3,L1);
F84 = F48;

% Surface 5
F56 = F13;
F65 = F31;
F58 = F18;
F85 = F81;

% Surface 6
F68 = F58;
F86 = F85;

F81 = F81/2;
F82 = F82/2;
F83 = F83/2;
F84 = F84/2;
F85 = F85/2;
F86 = F86/2;

% radiative coupling
R.R_12 = sigma_SB * A1 * epsilon_1 * epsilon_2 * F12;
R.R_1rad = 0;
R.R_13 = sigma_SB * A1 * epsilon_1 * epsilon_3 * F13;
R.R_14 = sigma_SB * A1 * epsilon_1 * epsilon_4 * F14;
R.R_15 = sigma_SB * A1 * epsilon_1 * epsilon_5 * F15;
R.R_16 = sigma_SB * A1 * epsilon_1 * epsilon_6 * F16;
R.R_18 = sigma_SB * A1 * epsilon_1 * epsilon_8 * F18;
R.R_23 = sigma_SB * A2 * epsilon_2 * epsilon_3 * F23;
R.R_24 = sigma_SB * A2 * epsilon_2 * epsilon_4 * F24;
R.R_25 = sigma_SB * A2 * epsilon_2 * epsilon_5 * F25;
R.R_26 = sigma_SB * A2 * epsilon_2 * epsilon_6 * F26;
R.R_28 = sigma_SB * A2 * epsilon_2 * epsilon_8 * F28;
R.R_rad3 = 0; % ?
R.R_rad4 = 0; % ?
R.R_rad2 = 0; % ?
R.R_rad6 = 0; % ?
R.R_rad8 = 0; % ?
R.R_rad5 = 0;
R.R_34 = sigma_SB * A3 * epsilon_3 * epsilon_4 * F34;
R.R_35 = sigma_SB * A3 * epsilon_3 * epsilon_5 * F35;
R.R_36 = sigma_SB * A3 * epsilon_3 * epsilon_6 * F36;
R.R_38 = sigma_SB * A3 * epsilon_3 * epsilon_8 * F38;
R.R_45 = sigma_SB * A4 * epsilon_4 * epsilon_5 * F45;
R.R_46 = sigma_SB * A4 * epsilon_4 * epsilon_6 * F46;
R.R_48 = sigma_SB * A4 * epsilon_4 * epsilon_8 * F48;
R.R_56 = sigma_SB * A5 * epsilon_5 * epsilon_6 * F56;
R.R_58 = sigma_SB * A5 * epsilon_5 * epsilon_8 * F58;
R.R_10 = sigma_SB * A1 * epsilon_MLI;
R.R_20 = sigma_SB * A2 * epsilon_MLI;
R.R_30 = sigma_SB * A3 * epsilon_MLI;
R.R_40 = sigma_SB * A4 * epsilon_MLI;
R.R_50 = sigma_SB * A5 * epsilon_MLI;
R.R_60 = sigma_SB * A6_ext * epsilon_MLI;
R.R_rad0 = sigma_SB * A_rad_tot * eps_rad;

% Add MLI on surface 3
R.R_3int3ext = sigma_SB * A3 * epsilon_MLI;
R.R_1int1ext = sigma_SB * A1 * epsilon_MLI;
R.R_2int2ext = sigma_SB * A2 * epsilon_MLI;
R.R_4int4ext = sigma_SB * A4 * epsilon_MLI;
R.R_5int5ext = sigma_SB * A5 * epsilon_MLI;
R.R_6int6ext = sigma_SB * A6_int * epsilon_MLI;


C.C_1int1ext = k_honeycomb*(A1/(l_str_tot));
C.C_2int2ext = k_honeycomb*(A2/(l_str_tot));
C.C_3int3ext = k_honeycomb*(A3/(l_str_tot));
C.C_4int4ext = k_honeycomb*(A4/(l_str_tot));
C.C_5int5ext = k_honeycomb*(A5/(l_str_tot));
C.C_6int6ext = k_honeycomb*(A6_ext/(l_str_tot));

% Conductive Coupling
% Conductive Coupling
nc = 100;
C1 = k_str*l_str*L1/(L2/2);
C2 = k_str*l_str*L1/(L3/2);
C_cont = nc*A1;
C.C_12 = (1/C1 + 1/C2 + 1/C_cont)^(-1);
C1 = k_str*l_str*L2/(L1/2);
C3 = k_str*l_str*L2/(L3/2);
C_cont = nc*A3;
C.C_13 = (1/C1 + 1/C3 + 1/C_cont)^(-1);
C1 = k_str*l_str*L1/(L2/2);
C4 = k_str*l_str*L1/(L3/2);
C_cont = nc*A4;
C.C_14 = (1/C1 + 1/C4 + 1/C_cont)^(-1);
C1 = k_str*l_str*L2/(L1/2);
C5 = k_str*l_str*L2/(L3/2);
C_cont = nc*A5;
C.C_15 = (1/C1 + 1/C5 + 1/C_cont)^(-1);
C.C_16 = 0;
C2 = k_str*l_str*L3/(L1/2);
C3 = k_str*l_str*L3/(L2/2);
C_cont = nc*A3;
C.C_23 = (1/C2 + 1/C3 + 1/C_cont)^(-1);
C2 = k_str*l_str*L3/(L1/2);
C5 = k_str*l_str*L3/(L2/2);
C_cont = nc*A5;
C.C_25 = (1/C2 + 1/C5 + 1/C_cont)^(-1);
C.C_24 = 0;
C.C_26 = C.C_12;
C.C_34 = C.C_23;
C.C_35 = 0;
C.C_36 = C.C_13;
C.C_36 = 0; %%%%%%%%%%%%%%%%%%%%
C.C_45 = C.C_34;
C.C_46 = C.C_14;
C.C_56 = C.C_13;

C.C_18 = (1/(k_str*l_str*L1*(+1/(L3/2))) + 1/(nc*l_str*L1))^(-1);
C.C_28 = 0;
C.C_38 = (1/(k_str*l_str*L3*(+1/(L1/2))) + 1/(nc*l_str*L3))^(-1);
C.C_48 = 0;
C.C_58 = (1/(k_str*l_str*L3*(+1/(L1/2))) + 1/(nc*l_str*L3))^(-1);
C.C_68 = C.C_18;
% to tune:
R_HP_min = 0.04;
R_HP_max = 0.76;
C_HP_max = 1/R_HP_min;
C_HP_min = 1/R_HP_max;
C_TS = 5;
C.C_1rad = 0;
C.C_2rad = 10;
C.C_3rad = 0;    % opposite surface
C.C_4rad = 0;
C.C_6rad = 10;
C.C_8rad = 10;
C.C_5rad = 3;


l_str_tot = 20e-3; % do sensitivity analysis
% C.C_6rad = k_honeycomb*(A_rad_tot/(l_str_tot));

% C.C_6rad = 2*C_TS;
% C.C_5rad = C.C_5rad + 4*C_HP_max;
% C.C_8rad = 6*C_TS;
% External fluxes 
q_Sun = q_sun_earth;
q_alb = q_Sun*F_earth*a_earth;
q_Earth = F_earth*sigma_SB*T_earth^4*epsilon_Earth;
q_alb_enc = q_sun_enc * F_enc_min * a_enc;
q_alb_sat = q_sun_enc * F_sat * a_sat;
% Angles with external fluxes

theta_6Sun = 24.7 * pi/180;
theta_3Sun = 90 - theta_6Sun;

Q_ext_hot = zeros(9,1);
% HOT CASE 1 : cameras towards Earth, face 5 sees the Sun
% HOT CASE 2: cameras towards Earth, antenna towards Sun
% HOT CASE 3: HGA (1) towards Earth, RTGs (6) towards Sun
% HOT CASE 4: FIRST FLYBY, PL towards Earth
% HOT CASE 5: FIRST FLYBY, HGA towards Earth
% HOT CASE 6: second FLYBY  PERICENTER, cameras towards Earth, radiators on face 6
hot_case = 6;

switch hot_case
    case 1
        theta_5Sun = 15 * pi/180;
        theta_6Sun = 75 * pi/180;
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
    case 4
        % first fly by: point cameras towards Earth
        theta_S1 = deg2rad(65.63);
        theta_S2 = deg2rad(70);
        theta_S3 = deg2rad(25);
        Q_ext_hot(1) = q_Sun * A1 * alpha_MLI* cos(theta_S1);
        Q_ext_hot(2) = 0;
        Q_ext_hot(3) = q_Earth * epsilon_MLI * A3 + q_alb * A3 * alpha_MLI + q_Sun * A3 * alpha_MLI* cos(theta_S3);
       case 5
        % first fly by: point HGA towards Earth
        theta_S5 = deg2rad(65.63);
        theta_S2 = deg2rad(70);
        theta_S1 = deg2rad(25);
        Q_ext_hot(1) = q_Sun * A1 * alpha_MLI* cos(theta_S1) + q_Earth * epsilon_MLI * A1 + q_alb * A1 * alpha_MLI;
        Q_ext_hot(2) = q_Sun * A2 * alpha_MLI* cos(theta_S2);
        Q_ext_hot(5) = q_Sun * A5_ext * alpha_MLI* cos(theta_S5);
        % Q_ext_hot(7) = q_Sun * A_rad_tot * alpha_louv_open* cos(theta_S5);
    case 6
        theta_S5 = deg2rad(9);
        theta_S1 = deg2rad(83.1);
        Q_ext_hot(1) = q_Sun * A1 * alpha_MLI* cos(theta_S1);
        Q_ext_hot(5) = q_Sun * A5 * alpha_MLI* cos(theta_S5);
        Q_ext_hot(3) = q_Earth * epsilon_MLI * A3 + q_alb * A3 * alpha_MLI;
end
% Initial condition
T0 = 293;

% solve
%%% Internal dissipation power
Q_diss_hot = zeros(8,1);
Q_diss_hot(5) =  Q_hot/2;
Q_diss_hot(2) =  Q_hot/2;

%% RTG
T_RTG = 170 + 273.15;
A_rtg_radiation = 580800e-6;
epsilon_RTG = 1;
% R.R_6RTG = sigma_SB * A_rtg_radiation * epsilon_MLI*epsilon_RTG;
R_6RTG = linspace(0,1e-8,40);
C_6RTG = linspace(0,5,50);


%%

% Initial condition
T0 = 293;

%%% SOLVE THE SYSTEM 
Clamped = 1;
T_guess = 273*ones(14,1);
T_land_hot_RTG = zeros(length(R_6RTG),length(C_6RTG),14);
for i = 1:length(R_6RTG)
    R.R_6RTG = R_6RTG(i);
    for k =1:length(C_6RTG)
        C.C_6RTG = C_6RTG(k);
options = optimoptions('fsolve','MaxFunctionEvaluations',50000,'Maxiterations',50000);
T_land_hot_RTG(i,k,:) = fsolve(@(T) HeatBalance_Lander_RTG(T, R, C, Q_ext_hot , Q_diss_hot, Clamped, T_RTG), T_guess, options);
    end
end

figure
surf(R_6RTG,C_6RTG,T_land_hot_RTG(:,:,8)'-273.15)
title('T8')

figure
surf(R_6RTG,C_6RTG,T_land_hot_RTG(:,:,2)'-273.15)
title('T2')

%% Cold case
% Power
P_budget_cold = 211.14;
P_input_TMTC_cold = 0;       % ask Antoine
P_diss_TMTC_cold =0;      % ask Antoine
% add batteries ...
We = 239;
Wt = 1589;
Q_cold = P_budget_cold-P_input_TMTC_cold+P_diss_TMTC_cold;

Q_diss_cold = zeros(8,1);
% Q_diss_cold(5) = Q_cold;

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
    theta_3Enc =deg2rad(50);
    theta_1Enc =deg2rad(40);
    theta_1Sat = 0; % CHANGE!
    Q_ext_cold(3) = q_Enc_orbit * A3 * epsilon_MLI*cos(theta_3Enc) ;
    Q_ext_cold(1) = q_Sat * A1 * epsilon_MLI*cos(theta_1Sat) + q_Enc_orbit * A1 * epsilon_MLI*cos(theta_1Sat) ;
    Clamped = 1;
    switch Clamped
        case 0
            % during landing: peak load. no eclipse !!
            C.C_2rad = 10;
            C.C_8rad = 12;
            theta_3Enc =deg2rad(10);
    theta_6Enc =deg2rad(80);
    theta_6Sat = 0; % CHANGE!
    theta_1Sun = deg2rad(1);
    Q_ext_cold(3) = q_Enc_orbit * A3 * epsilon_MLI*cos(theta_3Enc) + q_alb_enc * A3 * alpha_MLI*cos(theta_3Enc) ;
    Q_ext_cold(1) = q_sun_enc * A1 * alpha_MLI*cos(theta_1Sun) ;
    Q_ext_cold(6) = q_Enc_orbit * A6_ext * epsilon_MLI*cos(theta_6Enc) + q_alb_enc * A6_ext * alpha_MLI*cos(theta_6Enc) + ...
        q_Sat * A6_ext * epsilon_MLI*cos(theta_6Sat) + q_alb_sat * A6_ext * alpha_MLI*cos(theta_6Sat); 
Q_ext_cold(8) = q_Enc_orbit * A_rad_tot * eps_rad*cos(theta_6Enc) + q_alb_enc * A_rad_tot * alpha_louv_closed*cos(theta_6Enc) + ...
        q_Sat * A_rad_tot * eps_rad*cos(theta_6Sat) + q_alb_sat * A_rad_tot * alpha_louv_closed*cos(theta_6Sat); 
            % close louvers and compute again thermal couplings
            perc_open_rad = 3.5/n_rad; % percentage of open radiators
            eps_rad = perc_open_rad * eps_louv_open + (1-perc_open_rad) * eps_louv_closed;
            % eps_rad = eps_louv_closed;
            % radiative coupling
            R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
            C.C_2rad = 40;
            Q_diss_cold(6) = 344; % during landing --> check batteries..
        case 1
            C.C_6rad = 0;
            C.C_2rad = 3;
            C.C_5rad= 0; 
            C.C_8rad = 0;
            % clamped but in eclipse
            perc_open_rad = 0; % percentage of open radiators
            eps_rad = perc_open_rad * eps_louv_open + (1-perc_open_rad) * eps_louv_closed;
            % eps_rad = eps_louv_closed;
            % radiative coupling
            R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
            Q_shunt = 0;
            Q_diss_cold(2) = 100 + Q_shunt; % SAFE --> to change
    end
    case 2
        perc_open_rad = 2.1/n_rad; % percentage of open radiators
        eps_rad = perc_open_rad * eps_louv_open + (1-perc_open_rad) * eps_louv_closed;
        % eps_rad = eps_louv_closed;
        % radiative coupling
        Q_heaters = 0;
        C.C_2rad = 10;
        % Q_diss_cold(2) = 4;
        Q_diss_cold(8) = Q_heaters; % --> heaters
        Q_diss_cold(5) = Q_cold/2;
        Q_diss_cold(2) = Q_cold/2;
        R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
        theta_4Sat = deg2rad(20);
        theta_6Sat = deg2rad(20);
        theta_2Enc = 0;
        Clamped = 0;
        Q_ext_cold(2) = q_Enc_ground * A2 * epsilon_MLI * cos(theta_2Enc);
        Q_ext_cold(4) = q_Sat * A4*epsilon_MLI*cos(theta_4Sat);
        Q_ext_cold(1) = q_Sat * A1 *epsilon_MLI*cos(theta_6Sat);
end

%%% SOLVE THE SYSTEM 

%% SOLVE THE SYSTEM 
T_guess = 273*ones(14,1);
T_land_cold_RTG = zeros(length(R_6RTG),length(C_6RTG),14);
for i = 1:length(R_6RTG)
    R.R_6RTG = R_6RTG(i);
    for k =1:length(C_6RTG)
        C.C_6RTG = C_6RTG(k);
options = optimoptions('fsolve','MaxFunctionEvaluations',50000,'Maxiterations',50000);
T_land_cold_RTG(i,k,:) = fsolve(@(T) HeatBalance_Lander_RTG(T, R, C, Q_ext_cold , Q_diss_cold, Clamped, T_RTG), T_guess, options);
    end
end

figure
surf(R_6RTG,C_6RTG,T_land_cold_RTG(:,:,8)'-273.15)
title('T8')

figure
surf(R_6RTG,C_6RTG,T_land_cold_RTG(:,:,2)'-273.15)
title('T2')