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

% External fluxes 
q_Sun = q_sun_earth;
q_alb = q_Sun * F_earth * a_earth;
q_Earth = F_earth * sigma_SB * T_earth^4 * epsilon_Earth;

%%% Spacecraft %%%

% Area
L1 = 3.8;
L2 = 4.51;
L3 = 2.7;

% OLD = 
% L1 = 3.8;
% L2 = 4.51;
% L3 = 2.7;
nc = 100;


A1 = L1*L2;
A6_tot = A1;
A6_int = A6_tot;
A2 = L1*L3;
A4 = A2;
A3 = L2*L3;
A5_tot = A3;
A5_int = A5_tot;
D_ant = 3.8; % diameter of antenna [m]
A_ant = pi*(D_ant/2)^2; 
A1_ext = A1 - A_ant; % external area of surface 1 not covered by the antenna

%%% Structure properties (Al-5056-O)
k_str = 117;
% l_str = 0.02; %%%%%%%%%??
% l_str = 0.002; 
l_str = 0.002; % HYP : skin thickness of 5 mm

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
epsilon_15 = 0.874; 
% epsilon_5 = 0.034;
% epsilon_6 = 0.034;
% epsilon_15 = 0.034;
% epsilon_1 = 0.034;
% epsilon_2 = 0.034;
epsilon_3 = 0.034;
% epsilon_4 = 0.034;
% epsilon_MLI = 0.03;
epsilon_MLI = 0.02; %%%%%%%%%%%%%%%%%%%%%%
% alpha_MLI = 0.08 (Kapton, silvered, aluminum oxide coated, 1 mil)
%           = 0.11 Kapton, aluminized, silicon oxide coated, 1 mil
%           = 0.12 typical
alpha_MLI = 0.11; 

% Antenna painted white: Vita-var PV-100 white paint
% alpha_ant = 0.22; % from 0.17 to 0.27
% epsilon_ant = 0.82; % do sensitivity analysis up to 0.9
% antenna coating : PCBZ (book says 0.16 to 0.24 )
alpha_ant = 0.2; % BOL
% alpha_ant_EOL = 0.4;
epsilon_ant = 0.86;

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
A_rad_one = 1 *536*397e-6; % area one radiator.
% SENER LOUVER: 0.2128 m^2. consider also the structure... ok 0.2m^2
n_rad_5 = 15; %% can be changed
n_rad_6 = 10;
A_rad_tot_5 = A_rad_one*n_rad_5;% total area of radiators
A_rad_tot_6 = A_rad_one*n_rad_6;
A5_ext = A5_tot-A_rad_tot_5; % area on surface 5 not covered by radiators
A6_ext = A6_tot-A_rad_tot_6;
A_rad_tot = A_rad_tot_5 + A_rad_tot_6;
% RHU or heaters
n_RHU = 0; % you can change this!
P_RHU = 1; % 1 W is the thermal power generated by each RHU 


rho_rad = 6; % [kg/m^2] --> can be reduced to 1.5
m_louv = 1.07; % kg
mass_rad_5 = rho_rad*A_rad_tot_5;
mass_louv_5 = + m_louv * n_rad_5;
mass_rad_6 = rho_rad*A_rad_tot_6;
mass_louv_6 = + m_louv * n_rad_6;
mass_louv_tot = mass_louv_5 + mass_louv_6 ;
mass_rad_tot = mass_rad_5 + mass_rad_6 ;

rho_MLI = 0.870; % [kg/m^2]
A_MLI = A1_ext + A2 + A3 + A4 + A5_ext + A6_ext;
mass_MLI = rho_MLI*A_MLI;
%% PAYLOADS
% Q required by heaters/ RHU in the cold case
% NAC
A_NAC = 612863e-6;
eps_NAC = 0.031; % aluminum (ASMAD)
T_NAC= - 70 + 15 + 273.15;
Q_NAC_req = A_NAC*sigma_SB*T_NAC^4*eps_NAC;

% WAC
A_WAC = 983474e-6;
eps_WAC = 0.031; % aluminum (ASMAD)
T_WAC= - 30 + 15 + 273.15;
Q_WAC_req = A_WAC*sigma_SB*T_WAC^4*eps_WAC;

% TES
eps_TES = 0.031; % aluminum (ASMAD)
l1 = 130e-3; l2 = 180e-3; l3 = 180e-3;
A_TES = l1*l2 + 2*l2*l3 + 2*l1*l3;
T_TES = 10 + 15 + 273.15;
Q_TES_req = A_TES*sigma_SB*T_TES^4*eps_TES; % heat needed to warm up TES

% Laser altimeter
l1 = 600e-3;
l2 =  400e-3; l3 = 250e-3;
A_LA = l2*l3 + 2*l1*l3 + 2*l1*l2;
eps_LA = 0.031; % aluminum (ASMAD)
T_LA = - 20 + 15 + 273.15;
Q_LA_req = A_LA*sigma_SB*T_LA^4*eps_LA;

% Radar sounder
A_RS = 20341*2*(10^(-6))+21834e-6;
eps_RS = 0.031; % aluminum (ASMAD)

T_RS = - 40 + 15 + 273.15;
Q_RS_req = A_RS*sigma_SB*T_RS^4*eps_RS;

% sum
Q_pl_tot = Q_NAC_req + Q_WAC_req + Q_TES_req + Q_RS_req + Q_LA_req;
fprintf("Total power required by PL in the cold case [W]: %f ",Q_pl_tot)

% Compute area radiators for pl (hyp: power dissipated by electronics, not
% by pl)
alpha_al = 0.16; % book
% NAC
A_NAC_hot = pi*(390e-3/2)^2 -pi*(140e-3/2)^2 ;
Q_NAC = 7; % power budget
Q_hot_NAC = q_Earth * eps_NAC*A_NAC_hot + q_alb*A_NAC_hot*alpha_al + Q_NAC;
T_NAC_hot = 40 - 15 + 273.15;
Q_e_NAC =  A_NAC*sigma_SB*T_NAC_hot^4*eps_NAC;
Q_rad_NAC = Q_hot_NAC - Q_e_NAC;
A_rad_NAC  = Q_rad_NAC/(sigma_SB*T_NAC_hot^4*eps_rad);

% WAC
A_WAC_hot = pi*(760e-3/2)^2 -pi*( 600e-3/2)^2 ;
Q_WAC = 4;
Q_hot_WAC = q_Earth * eps_WAC*A_WAC_hot + q_alb*A_WAC_hot*alpha_al + Q_WAC;
T_WAC_hot = 40 - 15 + 273.15;
Q_e_WAC =  A_WAC*sigma_SB*T_WAC_hot^4*eps_WAC;
Q_rad_WAC = Q_hot_WAC - Q_e_WAC;
A_rad_WAC  = Q_rad_WAC/(sigma_SB*T_WAC_hot^4*eps_rad);

% TES
A_TES_hot = 23400E-6;
Q_TES = 18;
Q_hot_TES = q_Earth * eps_TES*A_TES_hot + q_alb*A_TES_hot*alpha_al + Q_TES;
T_TES_hot = 40 - 15 + 273.15;
Q_e_TES =  A_TES*sigma_SB*T_TES_hot^4*eps_TES;
Q_rad_TES = Q_hot_TES - Q_e_TES;
A_rad_TES  = Q_rad_TES/(sigma_SB*T_TES_hot^4*eps_rad);

% LASER ALTIMETER
A_LA_hot = 600e-3*250e-3;
Q_LA = 52;
Q_hot_LA = q_Earth * eps_TES*A_LA_hot + q_alb*A_LA_hot*alpha_al + Q_LA;
% T_LA_hot = 65 - 15 + 273.15;
T_LA_hot = 25 + 273.15;
Q_e_LA =  A_TES*sigma_SB*T_LA_hot^4*eps_LA;
Q_rad_LA = Q_hot_LA - Q_e_LA;
A_rad_LA  = Q_rad_LA/(sigma_SB*T_LA_hot^4*eps_rad);

A_rad_pl_tot = A_rad_NAC + A_rad_WAC + A_rad_TES + A_rad_LA; % AREA ON NODE 6 
fprintf("Total area radiator required by PL in the hot case [m^2]: %f ",A_rad_pl_tot)
Q_pl_budget = Q_NAC + Q_WAC + Q_TES + Q_LA;
A6_ext = A6_ext - A_rad_pl_tot;
%%
% Power
P_budget_hot = 347.56;
P_input_TMTC_h = 29.93; % ask Antoine %%%%%%%%%%%%%
P_diss_TMTC_h = 23.96;  % ask Antoine %%%%%%%%%%%%%
% add batteries ...
We = 456;              % electrical power from RTG
Wt = 3040;             % thermal power from RTG
Q_RHU_added = 0;
Q_pl_budget_mode = 74;
Q_hot = P_budget_hot - P_input_TMTC_h+P_diss_TMTC_h - Q_pl_budget_mode + Q_RHU_added; 
Q_hot = P_budget_hot - P_input_TMTC_h+P_diss_TMTC_h - Q_pl_budget_mode; 
Q_shunt = We - Q_hot; % check if a shunt can dissipate this power

% view factor
% % Surface 1
% F12 = VF_PerpRec(L3,L2,L1);
% F21 = F12*A1/A2;
% F13 = VF_PerpRec(L3,L1,L3);
% F31 = F13*A1/A3;
% F14 = VF_PerpRec(L3,L2,L1);
% F41 = VF_PerpRec(L2,L3,L1);
% F15_tot = F13;
% F51_tot = F15_tot*A5_tot/A1;
% F1rad = F15_tot*A_rad_tot/A5_tot;
% Frad1 = F1rad*A1/A_rad_tot; 
% F15 = F15_tot*A5/A5_tot;
% F51 = F15*A1/A5;
% F16 = VF_ParallelEqualRec(L3,L1,L2);
% F61 = F16;
H15 = 2600e-3;

F12 = 1*H15/L2;
F21 = VF_PerpRec(L2,L3,L1);
F13 = VF_PerpRec(L3,L1,L2-H15);
F31 = F13*(L2-H15)/L3;
F14 = VF_PerpRec(L3,L2-H15,L1);
F41 = F14*(L2-H15)/L3;
F15 = F13;
F51 = F31;
F16 = VF_ParallelEqualRec(L3,L1,L2-H15);
F61 = F16;
F115 = VF_PerpRec(L3,L2-H15,L1);
% Surface 2 
% F23 = VF_PerpRec(L2,L3,L1);
% F32 = F23*A2/A3;
% F24 = VF_ParallelEqualRec(L2,L1,L3);
% F42 = F24*A2/A4;
% F25_tot = VF_PerpRec(L2,L3,L1);
% F52_tot = F25_tot*A5_tot/A2;
% F2rad = F25_tot*A_rad_tot/A5_tot;
% Frad2 = F2rad*A2/A_rad_tot; 
% F25 = F25_tot*A5/A5_tot;
% F52 = F25*A2/A5;
% F26 = VF_PerpRec(L2,L3,L1);
% F62 = F26*A2/A6;
F32 = 1*H15/L2;
F23 = VF_PerpRec(L2,L1,L3);
F215 = VF_ParallelEqualRec(H15,L1,L3);
F315 = VF_PerpRec(L1,L2-H15,L3);
F42 = 0;
F415 = VF_ParallelEqualRec(L2-H15,L1,L3);
F515 = F315;
F615 = F115;
F24 = 0;
F62 = F12;
F26 = F21;
F52 = F32;
Frad2 = 0;
F25 = F23;
F2rad = 0;

% Surface 3
% F34 = F32;
% F43 = F23;
% F35_tot = VF_ParallelEqualRec(L1,L2,L3);
% F53_tot = F35_tot*A5_tot/A3;
% F3rad = F35_tot*A_rad_tot/A5_tot;
% Frad3 = F3rad*A3/A_rad_tot; 
% F35 = F35_tot*A5/A5_tot;
% F53 = F35*A3/A5;
% F36 = F31;
% F63 = F13;
F34 = VF_PerpRec(L1,L2-H15,L3);
F43 = F34*(L2-H15)/L1;
F35_tot =  VF_ParallelEqualRec(L1,H15,L3);
F53_tot = F35_tot;
F3rad = 0;
Frad3 = 0; 
F35 =  VF_ParallelEqualRec(L1,H15,L3);
F53 = F35_tot;
F36 = F31;
F63 = F13;

% Surface 4
% F45_tot = F43;
% F54_tot = F45_tot*A5_tot/A4;
% F4rad = F45_tot*A_rad_tot/A5_tot;
% Frad4 = F4rad*A4/A_rad_tot; 
% F45 = F45_tot*A5/A5_tot;
% F54 = F45*A4/A5;
% F46 = F41; 
% F64 = F14;

F45_tot = 0;
F54_tot = 0;
F4rad = 0;
Frad4 = 0; 
F45 = F43;
F54 = F34;
F46 = F41;
F64 = F14;

% Surface 6
% F65_tot = F63;
% F56_tot = F65_tot*A5_tot/A6;
% F6rad = F65_tot*A_rad_tot/A5_tot;
% Frad6 = F6rad*A6/A_rad_tot; 
% F65 = F65_tot*A5/A5_tot;
% F56 = F65*A6/A5;

F65_tot = 0;
F56_tot = 0;
F6rad = 0;
Frad6 = 0; 
F65 = 0;
F56 = 0;
F1rad = 0;
% radiative coupling
R.R_12 = sigma_SB * A1 * epsilon_1*epsilon_2 * F12;
R.R_1rad = 0;
R.R_13 = sigma_SB * A1 * epsilon_1*epsilon_3 * F13;
R.R_14 = sigma_SB * A1 * epsilon_1*epsilon_4 * F14;
R.R_15 = sigma_SB * A1 * epsilon_1*epsilon_5 * F15;
R.R_16 = sigma_SB * A1 * epsilon_1*epsilon_6 * F16;
R.R_23 = sigma_SB * A2 * epsilon_2*epsilon_3 * F23;
% R.R_24 = sigma_SB * A2 * epsilon_2*epsilon_4 * F24;
R.R_24 = 0;
R.R_25 = sigma_SB * A2 * epsilon_2*epsilon_5 * F25;
R.R_26 = sigma_SB * A2 * epsilon_2*epsilon_6 * F26;
R.R_rad3 = 0; % ?
R.R_rad4 = 0; % ?
R.R_rad2 = 0; % ?
R.R_rad6 = sigma_SB * A_rad_tot_6 * epsilon_MLI; % ?
R.R_rad5 = sigma_SB * A_rad_tot_5 * epsilon_MLI;
R.R_34 = sigma_SB * A3 * epsilon_3*epsilon_4 * F34;
R.R_35 = sigma_SB * A3 * epsilon_3*epsilon_5 * F35;
R.R_36 = sigma_SB * A3 * epsilon_3*epsilon_6 * F36;
R.R_45 = sigma_SB * A4 * epsilon_4*epsilon_5 * F45;
R.R_46 = sigma_SB * A4 * epsilon_4*epsilon_6 * F46;
R.R_56 = sigma_SB * A5_int * epsilon_5*epsilon_6 * F56;

R.R_115 = sigma_SB * A1 * epsilon_1*epsilon_15 * F115;
R.R_215 = sigma_SB * A2 * epsilon_2*epsilon_15 * F215;
R.R_315= sigma_SB * A3 * epsilon_3*epsilon_15 * F315;
R.R_415= sigma_SB * A4 * epsilon_4*epsilon_15 * F415;
R.R_515= sigma_SB * A5_int * epsilon_5*epsilon_15 * F515;
R.R_615= sigma_SB * A6_int * epsilon_6*epsilon_15 * F615;

C.C_115 =( 1/(k_str*l_str*L1*(+1/(L3/2))) + 1/(nc*l_str*L1))^(-1);
C.C_215 = 0;
C.C_315 = (1/(k_str*l_str*L3*(+1/(L1/2)))+ 1/(nc*l_str*L3))^(-1);

C.C_315 = 0;
C.C_415  = 0;
C.C_515  = C.C_315 ;
C.C_615 = C.C_115;


R.R_1ant = sigma_SB * (A1-A1_ext) * epsilon_MLI; % ????? not sure about this. also conduction. and not only MLI
C_1ant_max = k_honeycomb* l_str_tot/(A1 - A1_ext);
C.C_1ant = C_1ant_max/10;  % hyp: diameter of contact antenna - structure is 1/10 of D antenna 
% C.C_1ant = 0;
R.R_1ant = sigma_SB * (A1-A1_ext) * epsilon_MLI; % ????? not sure about this. also conduction. and not only MLI
R.R_10 = sigma_SB * A1_ext * epsilon_MLI;
R.R_20 = sigma_SB * A2 * epsilon_MLI;
R.R_30 = sigma_SB * A3 * epsilon_MLI;
R.R_40 = sigma_SB * A4 * epsilon_MLI;
A4_clamped = A4 - 1.5*1.5;
R.R_40_clamped = sigma_SB * A4_clamped * epsilon_MLI;
R.R_50 = sigma_SB * A5_ext * epsilon_MLI;
R.R_60 = sigma_SB * A6_ext * epsilon_MLI;
% R.R_60 = sigma_SB * (A6) * epsilon_MLI;
R.R_ant0 = sigma_SB * A_ant * epsilon_ant;
R.R_rad0 = sigma_SB * A_rad_tot * eps_rad;
R.R_rad15 = 0;
% Add MLI on surface 3, 2, 6 
R.R_6int6ext = sigma_SB * A6_ext * epsilon_MLI;
R.R_3int3ext = sigma_SB * A3 * epsilon_MLI;
R.R_2int2ext = sigma_SB * A2 * epsilon_MLI;
R.R_1int1ext = sigma_SB * A1 * epsilon_MLI;
R.R_4int4ext = sigma_SB * A4 * epsilon_MLI;
R.R_5int5ext = sigma_SB * A5_ext * epsilon_MLI;

% Conduction normal to honeycomb panel
C.C_1int1ext = k_honeycomb*(A1/(l_str_tot));
C.C_2int2ext = k_honeycomb*(A2/(l_str_tot));
C.C_3int3ext = k_honeycomb*(A3/(l_str_tot));
C.C_4int4ext = k_honeycomb*(A4/(l_str_tot));
C.C_5int5ext = k_honeycomb*(A5_ext/(l_str_tot));
C.C_6int6ext = k_honeycomb*(A6_ext/(l_str_tot));

% Conductive Coupling
nc = 100;
C1 = k_str*l_str*L1/(L2/2);
C2 = k_str*l_str*L1/(L3/2);
C_cont = nc*l_str*L1;
C.C_12 = (1/C1 + 1/C2 + 1/C_cont)^(-1);
C1 = k_str*l_str*L2/(L1/2);
C3 = k_str*l_str*L2/(L3/2);
C_cont = nc*l_str*L2;
C.C_13 = (1/C1 + 1/C3 + 1/C_cont)^(-1);
C1 = k_str*l_str*L1/(L2/2);
C4 = k_str*l_str*L1/(L3/2);
C_cont = nc*l_str*L1;
C.C_14 = (1/C1 + 1/C4 + 1/C_cont)^(-1);
C1 = k_str*l_str*L2/(L1/2);
C5 = k_str*l_str*L2/(L3/2);
C_cont = nc*l_str*L2;
C.C_15 = (1/C1 + 1/C5 + 1/C_cont)^(-1);
C.C_16 = 0;
C2 = k_str*l_str*L3/(L1/2);
C3 = k_str*l_str*L3/(L2/2);
C_cont = nc*l_str*L3;
C.C_23 = (1/C2 + 1/C3 + 1/C_cont)^(-1);
C2 = k_str*l_str*L3/(L1/2);
C5 = k_str*l_str*L3/(L2/2);
C_cont = nc*l_str*L3;
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

l_str_tot = 20e-3;
C.C_5rad = k_honeycomb*(A_rad_tot_5/(l_str_tot));
C.C_6rad = k_honeycomb*(A_rad_tot_6/(l_str_tot));
C.C_15rad = 20;
% sensitivity analysis !!!!!!!!!!!!!!!!!!!!!!!!!!!
C.C_5rad = 500;
C.C_6rad = 40;
C.C_3rad = 100;

% to tune:
% heat pipes from Celsia heat pipes calculator
R_HP_min = 0.03;
R_HP_max = 0.76;
C_HP_max = 1/R_HP_min;
C_HP_min = 1/R_HP_max;
C_TS = 5; % conductive coupling thermal straps from HiPeR Flexlinks AIRBUS
C.C_1rad = 0;
C.C_2rad = 0;
% C.C_3rad = 0;    % to radiators on surface 6
C.C_4rad = 0;
% C.C_6rad = 0;
% C.C_5rad = 0; % HYP: HP both for electronics TMTC, batteries and OBDH

% Angles with external fluxes
theta_6Sun = 24.7 * pi/180;
theta_3Sun = 90 - theta_6Sun;

Q_ext_hot = zeros(8,1);
% HOT CASE 1 : cameras towards Earth, face 5 sees the Sun
% HOT CASE 2: cameras towards Earth, antenna towards Sun
% HOT CASE 3: HGA (1) towards Earth, RTGs (6) towards Sun
% HOT CASE 4: FIRST FLYBY, PL towards Earth
% HOT CASE 5: FIRST FLYBY, HGA towards Earth
hot_case = 4;

switch hot_case
    case 1
        theta_5Sun = deg2rad(20);
        theta_6Sun = pi/2 - theta_5Sun;
        Q_ext_hot(3) = q_Earth * epsilon_MLI*A3 + q_alb*A3*alpha_MLI ;
        Q_ext_hot(5) = q_Sun * A5 * alpha_MLI* cos(theta_5Sun); 
        Q_ext_hot(6) = q_Sun * A6 * alpha_MLI* cos(theta_6Sun);
    case 2
        theta_antSun = 0;
        Q_ext_hot(7) = q_Sun * A_ant * alpha_ant* cos(theta_antSun);
        Q_ext_hot(1) = q_Sun * A1_ext * alpha_MLI* cos(theta_antSun);
        Q_ext_hot(3) = q_Earth*epsilon_MLI*A3 + q_alb*A3*alpha_MLI ;
      
    case 3
        Q_ext_hot(6) = q_Sun * A6 * alpha_MLI* cos(theta_6Sun);
        Q_ext_hot(3) = q_Sun * A3 * alpha_MLI* cos(theta_3Sun);
        Q_ext_hot(1) = q_Earth * epsilon_MLI * A1_ext + q_alb * A1_ext * alpha_MLI ;
        Q_ext_hot(7) = q_Earth * epsilon_ant * A_ant + q_alb * A_ant * alpha_ant ;
    case 4
        % first fly by: point cameras towards Earth
        theta_S1 = deg2rad(65.63);
        theta_S2 = deg2rad(70);
        theta_S3 = deg2rad(25);
        Q_ext_hot(1) = q_Sun * A1_ext * alpha_MLI* cos(theta_S1);
        Q_ext_hot(2) = q_Sun * A2 * alpha_MLI* cos(theta_S2);
        Q_ext_hot(3) = q_Earth * epsilon_MLI * A3 + q_alb * A3 * alpha_MLI + q_Sun * A3 * alpha_MLI* cos(theta_S3);
        Q_ext_hot(7) = q_Sun * A_ant * alpha_ant* cos(theta_S1);
       case 5
        % first fly by: point HGA towards Earth
        theta_S5 = deg2rad(65.63);
        theta_S2 = deg2rad(70);
        theta_S1 = deg2rad(25);
        Q_ext_hot(1) = q_Sun * A1_ext * alpha_MLI* cos(theta_S1) + q_Earth * epsilon_MLI * A1_ext + q_alb * A1_ext * alpha_MLI;
        Q_ext_hot(2) = q_Sun * A2 * alpha_MLI* cos(theta_S2);
        Q_ext_hot(7) = q_Sun * A_ant * alpha_ant* cos(theta_S1) + q_Earth * epsilon_ant * A_ant + q_alb * A_ant * alpha_ant;
        Q_ext_hot(5) = q_Sun * A5_ext * alpha_MLI* cos(theta_S5);
        Q_ext_hot(8) = q_Sun * A_rad_tot_5 * alpha_louv_open* cos(theta_S5);
end

% Initial condition
T0 = 293;

% solve
%%% Internal dissipation power
Q_diss_hot = zeros(15,1);
Q_diss_hot(15) =  Q_hot/2;
Q_diss_hot(6) = Q_hot/2;
%%% SOLVE THE SYSTEM 
Clamped = 1;
T_guess = 273*ones(15,1);
options = optimoptions('fsolve','display','iter','MaxFunctionEvaluations',50000,'Maxiterations',50000);
T_orb_hot = fsolve(@(T) HeatBalance_Orbiter(T, R, C, Q_ext_hot , Q_diss_hot, Clamped), T_guess, options);

fprintf(['1 ',num2str(T_orb_hot(1)-273),' Celsius\n'])
fprintf(['2 ',num2str(T_orb_hot(2)-273),' Celsius\n'])
fprintf(['3 ',num2str(T_orb_hot(3)-273),' Celsius\n'])
fprintf(['4 ',num2str(T_orb_hot(4)-273),' Celsius\n'])
fprintf(['5 ',num2str(T_orb_hot(5)-273),' Celsius\n'])
fprintf(['6 ',num2str(T_orb_hot(6)-273),' Celsius\n'])
fprintf(['ant ',num2str(T_orb_hot(7)-273),' Celsius\n'])
fprintf(['rad ',num2str(T_orb_hot(8)-273),' Celsius\n'])
fprintf(['15',num2str(T_orb_hot(15)-273),' Celsius\n'])

% add mass and specific heat for transient


%% Cold case
% Power
P_budget_cold = 258; % RW desaturation
P_input_TMTC_cold = 0;
P_diss_TMTC_cold = 0;
% add batteries ...
We = 389;
We_TCS = Q_pl_tot; % power required by pl
We_av = We - We_TCS; % power available internally
Wt = 2596;
Q_cold = P_budget_cold-P_input_TMTC_cold+P_diss_TMTC_cold;

% close louvers and compute again thermal couplings
eps_rad = eps_louv_closed;

% radiative coupling

R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
C.C_6rad = 1.5; % ? sensitivity analysis
C.C_5rad = 0; % ? sensitivity analysis
% C.C_13 = C.C_13 + 30;
C.C_3rad = 5;
C.C_3rad = 3.8;
C.C_15rad = 2;
% External fluxes
% IR Heat fluxes for Saturn and Enceladus
q_Sat = F_sat*sigma_SB*T_Sat^4*epsilon_sat;
q_Enc = F_enc_min * sigma_SB * T_Enc^4 *epsilon_Enc;

% HYP: nadir pointing with cameras
% face 3 with cameras points Enceladus and receive IR from Enceladus
% Saturn ? HYP: face 4 (CHANGE!)
P_added_6 = 150; % from RTG (look Curiosity, from pumped fluid loop)
P_VRHU_6 = 0; % for batteries
P_VRHU_3 = 0; % for internal PL
% Q_diss_cold = zeros(8,1);
% Q_diss_cold(6) = P_added_6  + We_av*perc + P_VRHU_6;
% Q_diss_cold(3) = P_VRHU_3;
% Q_diss_cold(5) = We_av*(1-perc);

% External fluxes
theta_3Enc = 0;
theta_6Sat = 0; % CHANGE!
Q_ext_cold = zeros(11,1);
Q_ext_cold(3) = q_Enc * A3 * epsilon_MLI*cos(theta_3Enc) ;
Q_ext_cold(6) = q_Sat * A6_ext * epsilon_MLI*cos(theta_6Sat) ;
Q_ext_cold(8) = q_Sat * A_rad_tot_6 * eps_rad*cos(theta_6Sat);

%%% SOLVE THE SYSTEM 
Clamped = 1;
if Clamped == 0
% P_added_6 = 150; % from RTG (look Cassini Reference)
P_added_6 = 230;
P_VRHU_6 = 0; % for batteries
P_VRHU_3 = 0; % for internal PL
perc = 0.7;
Q_diss_cold = zeros(15,1);
Q_diss_cold(6) = P_added_6 + We_av*perc + P_VRHU_6;
Q_diss_cold(3) = P_VRHU_3;
Q_diss_cold(15) = We_av*(1-perc);
else 
    if Clamped == 1
        P_added_6 = 225; % from RTG (look Cassini Reference)
        P_VRHU_6 = 0; % for batteries
        P_VRHU_3 = 0; % for internal PL
        perc = 0;
        Q_diss_cold = zeros(15,1);
        Q_diss_cold(6) = P_added_6 + We_av*perc + P_VRHU_6;
        Q_diss_cold(3) = P_VRHU_3;
        Q_diss_cold(15) = We_av*(1-perc);
    end
end
theta_Sun_enc = deg2rad(1); % Ask Antoine: min angle Sun-Earth during communications
Q_ext_cold(7) = q_sun_enc*A_ant*alpha_ant*cos(theta_Sun_enc);

T_guess = 273*ones(15,1);
options = optimoptions('fsolve','display','iter','MaxFunctionEvaluations',50000,'Maxiterations',50000);
T_orb_cold = fsolve(@(T) HeatBalance_Orbiter(T, R, C, Q_ext_cold , Q_diss_cold, Clamped), T_guess, options);

fprintf(['1 ',num2str(T_orb_cold(1)-273),' Celsius\n'])
fprintf(['2 ',num2str(T_orb_cold(2)-273),' Celsius\n'])
fprintf(['3 ',num2str(T_orb_cold(3)-273),' Celsius\n'])
fprintf(['4 ',num2str(T_orb_cold(4)-273),' Celsius\n'])
fprintf(['5 ',num2str(T_orb_cold(5)-273),' Celsius\n'])
fprintf(['6 ',num2str(T_orb_cold(6)-273),' Celsius\n'])
fprintf(['ant ',num2str(T_orb_cold(7)-273),' Celsius\n'])
fprintf(['rad ',num2str(T_orb_cold(8)-273),' Celsius\n'])
fprintf(['15 ',num2str(T_orb_cold(15)-273),' Celsius\n'])

% Conduction between surfaces
% To do:
% 2) change internal power
% 3) surface in contact with lander
% 4) External p/l ?
% sensitivity analysis, check properties
% T of RTG

%% Sensitivity analysis:

%% k honeycomb


k_honeycomb_vec = linspace(0.054,2.036,20);
T_orb_cold_vec2 = zeros(20,15);
T_orb_hot_vec2 = zeros(20,15);
for m = 1:20
    k_honeycomb = k_honeycomb_vec(m);
C_1ant_max = k_honeycomb* l_str_tot/(A1 - A1_ext);
% Conduction normal to honeycomb panel
C.C_1int1ext = k_honeycomb*(A1/(l_str_tot));
C.C_2int2ext = k_honeycomb*(A2/(l_str_tot));
C.C_3int3ext = k_honeycomb*(A3/(l_str_tot));
C.C_4int4ext = k_honeycomb*(A4/(l_str_tot));
C.C_5int5ext = k_honeycomb*(A5_ext/(l_str_tot));
C.C_6int6ext = k_honeycomb*(A6_ext/(l_str_tot));
T_guess = 273*ones(15,1);
options = optimoptions('fsolve','MaxFunctionEvaluations',50000,'Maxiterations',50000);
C.C_6rad = 1.5; % ? sensitivity analysis
C.C_5rad = 0; % ? sensitivity analysis
% C.C_13 = C.C_13 + 30;
C.C_3rad = 5;
C.C_15rad = 2;
eps_rad = eps_louv_closed;
R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
T_orb_cold_vec2(m,:) = fsolve(@(T) HeatBalance_Orbiter(T, R, C, Q_ext_cold , Q_diss_cold, Clamped), T_guess, options);
C.C_5rad = k_honeycomb*(A_rad_tot_5/(l_str_tot));
C.C_6rad = k_honeycomb*(A_rad_tot_6/(l_str_tot));
C.C_15rad = 20;
% sensitivity analysis !!!!!!!!!!!!!!!!!!!!!!!!!!!
C.C_5rad = 500;
C.C_6rad = 40;
C.C_3rad = 100;
eps_rad = eps_louv_open;
R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
T_orb_hot_vec2(m,:) = fsolve(@(T) HeatBalance_Orbiter(T, R, C, Q_ext_hot , Q_diss_hot, Clamped), T_guess, options);

end
figure
plot(k_honeycomb_vec,T_orb_cold_vec2(:,3)-273.15)
figure
plot(k_honeycomb_vec,T_orb_hot_vec2(:,3)-273.15)


figure
plot(k_honeycomb_vec,T_orb_cold_vec2(:,15)-273.15)
figure
plot(k_honeycomb_vec,T_orb_hot_vec2(:,15)-273.15)
%%
% conduction
k_honeycomb = 2.036;
C_1ant_max = k_honeycomb* l_str_tot/(A1 - A1_ext);
% Conduction normal to honeycomb panel
C.C_1int1ext = k_honeycomb*(A1/(l_str_tot));
C.C_2int2ext = k_honeycomb*(A2/(l_str_tot));
C.C_3int3ext = k_honeycomb*(A3/(l_str_tot));
C.C_4int4ext = k_honeycomb*(A4/(l_str_tot));
C.C_5int5ext = k_honeycomb*(A5_ext/(l_str_tot));
C.C_6int6ext = k_honeycomb*(A6_ext/(l_str_tot));

k_str_vec = linspace(88,229,20);
l_str_vec = linspace(0.005,0.004,20);
T_orb_cold_vec = zeros(20,20,15);
T_orb_hot_vec = zeros(20,20,15);
for k = 1:length(k_str_vec)
    for j = 1:length(l_str_vec)
        k_str =k_str_vec(k);
        l_str = l_str_vec(j);
        C.C_115 =( 1/(k_str*l_str*L1*(+1/(L3/2))) + 1/(nc*l_str*L1))^(-1);
        C.C_215 = 0;
        C.C_315 = (1/(k_str*l_str*L3*(+1/(L1/2)))+ 1/(nc*l_str*L3))^(-1);
        nc = 100;
C1 = k_str*l_str*L1/(L2/2);
C2 = k_str*l_str*L1/(L3/2);
C_cont = nc*l_str*L1;
C.C_12 = (1/C1 + 1/C2 + 1/C_cont)^(-1);
C1 = k_str*l_str*L2/(L1/2);
C3 = k_str*l_str*L2/(L3/2);
C_cont = nc*l_str*L2;
C.C_13 = (1/C1 + 1/C3 + 1/C_cont)^(-1);
C1 = k_str*l_str*L1/(L2/2);
C4 = k_str*l_str*L1/(L3/2);
C_cont = nc*l_str*L1;
C.C_14 = (1/C1 + 1/C4 + 1/C_cont)^(-1);
C1 = k_str*l_str*L2/(L1/2);
C5 = k_str*l_str*L2/(L3/2);
C_cont = nc*l_str*L2;
C.C_15 = (1/C1 + 1/C5 + 1/C_cont)^(-1);
C.C_16 = 0;
C2 = k_str*l_str*L3/(L1/2);
C3 = k_str*l_str*L3/(L2/2);
C_cont = nc*l_str*L3;
C.C_23 = (1/C2 + 1/C3 + 1/C_cont)^(-1);
C2 = k_str*l_str*L3/(L1/2);
C5 = k_str*l_str*L3/(L2/2);
C_cont = nc*l_str*L3;
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
T_guess = 273*ones(15,1);
options = optimoptions('fsolve','MaxFunctionEvaluations',50000,'Maxiterations',50000);
C.C_6rad = 1.5; % ? sensitivity analysis
C.C_5rad = 0; % ? sensitivity analysis
% C.C_13 = C.C_13 + 30;
C.C_3rad = 5;
C.C_15rad = 2;
eps_rad = eps_louv_closed;
R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
T_orb_cold_vec(k,j,:) = fsolve(@(T) HeatBalance_Orbiter(T, R, C, Q_ext_cold , Q_diss_cold, Clamped), T_guess, options);
C.C_5rad = k_honeycomb*(A_rad_tot_5/(l_str_tot));
C.C_6rad = k_honeycomb*(A_rad_tot_6/(l_str_tot));
C.C_15rad = 20;
% sensitivity analysis !!!!!!!!!!!!!!!!!!!!!!!!!!!
C.C_5rad = 500;
C.C_6rad = 40;
C.C_3rad = 100;
eps_rad = eps_louv_open;
R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
T_orb_hot_vec(k,j,:) = fsolve(@(T) HeatBalance_Orbiter(T, R, C, Q_ext_hot , Q_diss_hot, Clamped), T_guess, options);
    end
end

figure
surf(k_str_vec,l_str_vec,T_orb_cold_vec(:,:,1)-273.15)
figure
surf(k_str_vec,l_str_vec,T_orb_hot_vec(:,:,1)-273.15)
figure
surf(k_str_vec,l_str_vec,T_orb_cold_vec(:,:,2)-273.15)
figure
surf(k_str_vec,l_str_vec,T_orb_hot_vec(:,:,2)-273.15)

figure
surf(k_str_vec,l_str_vec,T_orb_cold_vec(:,:,3)-273.15)
figure
surf(k_str_vec,l_str_vec,T_orb_hot_vec(:,:,3)-273.15)


figure
surf(k_str_vec,l_str_vec,T_orb_cold_vec(:,:,15)-273.15)
figure
surf(k_str_vec,l_str_vec,T_orb_hot_vec(:,:,15)-273.15)

%%
% sensitivity analysis on alpha MLI and epsilon MLI --> find area radiators
theta_S1 = deg2rad(65.63);
theta_S2 = deg2rad(70);
theta_S3 = deg2rad(25);

% solve
%%% Internal dissipation power
Q_diss_hot = zeros(15,1);
Q_diss_hot(15) =  Q_hot/2;
Q_diss_hot(6) = Q_hot/2;
%%% SOLVE THE SYSTEM 
Clamped = 1;
y_guess = [273*ones(2,1);C.C_3rad;273*ones(12,1); A_rad_tot];
epsilon_MLI_vec =  linspace(0.005,0.05,20);
% epsilon_MLI_vec =  ones(3,1)*epsilon_MLI;
alpha_MLI_vec =  linspace(0.06,0.13,20);
Arad = zeros(length(epsilon_MLI_vec),length(alpha_MLI_vec),1);
T3 = zeros(length(epsilon_MLI_vec),length(alpha_MLI_vec),1);
% for k = 1:length(epsilon_MLI_vec)
%     eps_MLI = epsilon_MLI_vec(k);
%     R.R_rad6 = sigma_SB * A_rad_tot_6 * eps_MLI; % ?
%     R.R_rad5 = sigma_SB * A_rad_tot_5 * eps_MLI;
%     R.R_1ant = sigma_SB * (A1-A1_ext) * eps_MLI; % ????? not sure about this. also conduction. and not only MLI
%     R.R_10 = sigma_SB * A1_ext * eps_MLI;
%     R.R_20 = sigma_SB * A2 * eps_MLI;
%     R.R_30 = sigma_SB * A3 * eps_MLI;
%     R.R_40 = sigma_SB * A4 * eps_MLI;
%     A4_clamped = A4 - 1.5*1.5;
%     R.R_40_clamped = sigma_SB * A4_clamped * eps_MLI;
%     R.R_50 = sigma_SB * A5_ext * eps_MLI;
%     R.R_60 = sigma_SB * A6_ext * eps_MLI;
%     % Add MLI on surface 3, 2, 6 
%     R.R_6int6ext = sigma_SB * A6_ext * eps_MLI;
%     R.R_3int3ext = sigma_SB * A3 * eps_MLI;
%     R.R_2int2ext = sigma_SB * A2 * eps_MLI;
%     R.R_1int1ext = sigma_SB * A1 * eps_MLI;
%     R.R_4int4ext = sigma_SB * A4 * eps_MLI;
%     R.R_5int5ext = sigma_SB * A5_ext * eps_MLI;
%     for j = 1:20
%         alp_MLI = alpha_MLI_vec(j);
%         Q_ext_hot(1) = q_Sun * A1_ext * alp_MLI* cos(theta_S1);
%         Q_ext_hot(2) = q_Sun * A2 * alp_MLI* cos(theta_S2);
%         Q_ext_hot(3) = q_Earth * eps_MLI * A3 + q_alb * A3 * alp_MLI + q_Sun * A3 * alp_MLI* cos(theta_S3);
%         Q_ext_hot(7) = q_Sun * A_ant * alpha_ant* cos(theta_S1);
% 
%         options = optimoptions('fsolve','MaxFunctionEvaluations',50000,'Maxiterations',50000);
%         y_orb_hot = fsolve(@(y) OrbiterAreaRad(y, R, C, Q_ext_hot , Q_diss_hot, Clamped, A6_tot, sigma_SB, eps_MLI, eps_rad), T_guess, options);
%         Arad(j,k) = y_orb_hot(15);
%         T3(j,k) = y_orb_hot(3);
%     end
% end

% figure
% surf(epsilon_MLI_vec,alpha_MLI_vec, Arad)
% xlabel('epsilon')
% ylabel('alpha')
% zlabel('Arad')
% figure
% surf(epsilon_MLI_vec,alpha_MLI_vec,T3-273.15)
% xlabel('epsilon')
% ylabel('alpha')
% zlabel('T3')

% epsilon
Arad = zeros(length(epsilon_MLI_vec),1);
C3rad = zeros(length(epsilon_MLI_vec),1);
eps_rad = eps_louv_open;
R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
for k = 1:length(epsilon_MLI_vec)
    eps_MLI = epsilon_MLI_vec(k);
    R.R_rad6 = sigma_SB * A_rad_tot_6 * eps_MLI; % ?
    R.R_rad5 = sigma_SB * A_rad_tot_5 * eps_MLI;
    R.R_1ant = sigma_SB * (A1-A1_ext) * eps_MLI; % ????? not sure about this. also conduction. and not only MLI
    R.R_10 = sigma_SB * A1_ext * eps_MLI;
    R.R_20 = sigma_SB * A2 * eps_MLI;
    R.R_30 = sigma_SB * A3 * eps_MLI;
    R.R_40 = sigma_SB * A4 * eps_MLI;
    A4_clamped = A4 - 1.5*1.5;
    R.R_40_clamped = sigma_SB * A4_clamped * eps_MLI;
    R.R_50 = sigma_SB * A5_ext * eps_MLI;
    R.R_60 = sigma_SB * A6_ext * eps_MLI;
    % Add MLI on surface 3, 2, 6 
    R.R_6int6ext = sigma_SB * A6_ext * eps_MLI;
    R.R_3int3ext = sigma_SB * A3 * eps_MLI;
    R.R_2int2ext = sigma_SB * A2 * eps_MLI;
    R.R_1int1ext = sigma_SB * A1 * eps_MLI;
    R.R_4int4ext = sigma_SB * A4 * eps_MLI;
    R.R_5int5ext = sigma_SB * A5_ext * eps_MLI;
   
        alp_MLI = alpha_MLI;
        Q_ext_hot(1) = q_Sun * A1_ext * alp_MLI* cos(theta_S1);
        Q_ext_hot(2) = q_Sun * A2 * alp_MLI* cos(theta_S2);
        Q_ext_hot(3) = q_Earth * eps_MLI * A3 + q_alb * A3 * alp_MLI + q_Sun * A3 * alp_MLI* cos(theta_S3);
        Q_ext_hot(7) = q_Sun * A_ant * alpha_ant* cos(theta_S1);

        options = optimoptions('fsolve','MaxFunctionEvaluations',50000,'Maxiterations',50000);
        y_orb_hot = fsolve(@(y) OrbiterAreaRad(y, R, C, Q_ext_hot , Q_diss_hot, Clamped, A6_tot, sigma_SB, eps_MLI, eps_rad), y_guess, options);
        Arad(k) = y_orb_hot(15);
        C3rad(k) = y_orb_hot(3);
end

figure
plot(epsilon_MLI_vec, Arad)
xlabel('epsilon')
ylabel('Arad')
figure
plot(epsilon_MLI_vec,C3rad)
xlabel('epsilon')
ylabel('C3rad')

%% alpha
Arad = zeros(length(epsilon_MLI_vec),1);
C3rad = zeros(length(epsilon_MLI_vec),1);
eps_MLI = epsilon_MLI;
eps_rad = eps_louv_open;
R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
    for j = 1:20
        alp_MLI = alpha_MLI_vec(j);
        Q_ext_hot(1) = q_Sun * A1_ext * alp_MLI* cos(theta_S1);
        Q_ext_hot(2) = q_Sun * A2 * alp_MLI* cos(theta_S2);
        Q_ext_hot(3) = q_Earth * eps_MLI * A3 + q_alb * A3 * alp_MLI + q_Sun * A3 * alp_MLI* cos(theta_S3);
        Q_ext_hot(7) = q_Sun * A_ant * alpha_ant* cos(theta_S1);

        options = optimoptions('fsolve','MaxFunctionEvaluations',50000,'Maxiterations',50000);
        y_orb_hot = fsolve(@(y) OrbiterAreaRad(y, R, C, Q_ext_hot , Q_diss_hot, Clamped, A6_tot, sigma_SB, eps_MLI, eps_rad), y_guess, options);
        Arad(j) = y_orb_hot(15);
        C3rad(j) = y_orb_hot(3);
    end


figure
plot(alpha_MLI_vec,Arad)
xlabel('alpha')
ylabel('Arad')
figure
plot(alpha_MLI_vec,C3rad)

xlabel('alpha')
ylabel('C3rad')


%% Cold case --> change epsilon and see how much heat you need
y_guess = [273*ones(2,1);C.C_3rad;273*ones(12,1); P_added_6];
% epsilon
Wadded = zeros(length(epsilon_MLI_vec),1);
C3rad = zeros(length(epsilon_MLI_vec),1);
Q_ext_cold = zeros(11,1);
for k = 1:length(epsilon_MLI_vec)
    eps_MLI = epsilon_MLI_vec(k);
    R.R_rad6 = sigma_SB * A_rad_tot_6 * eps_MLI; % ?
    R.R_rad5 = sigma_SB * A_rad_tot_5 * eps_MLI;
    R.R_1ant = sigma_SB * (A1-A1_ext) * eps_MLI; % ????? not sure about this. also conduction. and not only MLI
    R.R_10 = sigma_SB * A1_ext * eps_MLI;
    R.R_20 = sigma_SB * A2 * eps_MLI;
    R.R_30 = sigma_SB * A3 * eps_MLI;
    R.R_40 = sigma_SB * A4 * eps_MLI;
    A4_clamped = A4 - 1.5*1.5;
    R.R_40_clamped = sigma_SB * A4_clamped * eps_MLI;
    R.R_50 = sigma_SB * A5_ext * eps_MLI;
    R.R_60 = sigma_SB * A6_ext * eps_MLI;
    % Add MLI on surface 3, 2, 6 
    R.R_6int6ext = sigma_SB * A6_ext * eps_MLI;
    R.R_3int3ext = sigma_SB * A3 * eps_MLI;
    R.R_2int2ext = sigma_SB * A2 * eps_MLI;
    R.R_1int1ext = sigma_SB * A1 * eps_MLI;
    R.R_4int4ext = sigma_SB * A4 * eps_MLI;
    R.R_5int5ext = sigma_SB * A5_ext * eps_MLI;
   % close louvers and compute again thermal couplings
eps_rad = eps_louv_closed;

% radiative coupling

R.R_rad0 = sigma_SB*A_rad_tot * eps_rad;
C.C_6rad = 1.5; % ? sensitivity analysis
C.C_5rad = 0; % ? sensitivity analysis
% C.C_13 = C.C_13 + 30;
C.C_3rad = 5;
C.C_15rad = 2;


Q_ext_cold(3) = q_Enc * A3 * eps_MLI*cos(theta_3Enc) ;
Q_ext_cold(6) = q_Sat * A6_ext * eps_MLI*cos(theta_6Sat) ;
Q_ext_cold(8) = q_Sat * A_rad_tot_6 * eps_rad*cos(theta_6Sat);

%%% SOLVE THE SYSTEM 
Clamped = 1;
if Clamped == 0

P_added_6 = 0;
P_VRHU_6 = 0; % for batteries
P_VRHU_3 = 0; % for internal PL
perc = 0.7;
Q_diss_cold = zeros(15,1);
Q_diss_cold(6) = P_added_6 + We_av*perc + P_VRHU_6;
Q_diss_cold(3) = P_VRHU_3;
Q_diss_cold(15) = We_av*(1-perc);
else 
    if Clamped == 1
        P_added_6 = 0; % to size !

        perc = 0;
        Q_diss_cold = zeros(15,1);
        Q_diss_cold(6) = P_added_6 + We_av*perc + P_VRHU_6;

        Q_diss_cold(15) = We_av*(1-perc);
    end
end
       options = optimoptions('fsolve','MaxFunctionEvaluations',50000,'Maxiterations',50000);
        y_orb_hot = fsolve(@(y) OrbiterHeaters(y, R,C, Q_ext_cold , Q_diss_cold, Clamped, sigma_SB, epsilon_MLI, eps_rad, A6_ext), y_guess, options);
        Wadded(k) = y_orb_hot(15);
        C3rad(k) = y_orb_hot(3);
end

figure
plot(epsilon_MLI_vec, Wadded)
xlabel('epsilon')
ylabel('Wadded')
figure
plot(epsilon_MLI_vec,C3rad)
xlabel('epsilon')
ylabel('C3rad')

%%

% fallo anche per il lander
% sistema k_honeycomb
% alpha MLI
% epsilon MLI
% --> trova radiatori e calore da aggiungere


% RTG