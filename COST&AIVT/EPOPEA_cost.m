%% EPOPEA COST ESTIMATION

clear; clc; close all;


%% Hardware cost - CER relationships - New SMAD (pag 311)
% MASSES UPDATED ATV 03/01

% PAYLOAD & ROBOTICS
nac = 26000; wac = 26000;
tes = 25000; radar = 33000; altimeter = 36200;
micro = 20000; es = 17000; eif = es;
hrms = 93400; msla = 72700;
sample = 78000; seism = 12000;
context = 11600; active_sample = 78000;

PL_cost = nac+wac+tes+radar+altimeter+micro+es+hrms+msla+sample+seism+...
    context+active_sample; % [$], FY 2025

% If we consider two different units where also the payload are integrated
% and tested (eg: qualification unit and flight model)
PL_cost = PL_cost*2;       % [$], FY 2025

% SUBSYSTEMS cost [FY 2010]

% Sampling orbiter
MO_stm = 1331.4 + 114.4; % Structure + thermal control mass [kg]
MO_adcs = 51.9; % ADCS mass [kg]
MO_eps = 139.7; % EPS mass [kg]
MO_rcs = 1578400; % Volume of the tank only for secondary prop. [cm3]
MO_me = 4.9; t_me_O = 50299; % Main engine mass and firing time [kg],[s]
MO_tmtc = 145.1; % Communication system mass [kg]
% Sampling lander
ML_stm = 81+51.2; % Structure + thermal control mass [kg]
ML_adcs = 20; % ADCS mass [kg]
ML_eps = 116.2; % EPS mass [kg]
ML_rcs = 8100; % Volume of the tank only for secondary prop. [cm3]
ML_tmtc = 32.3; % Communication system mass [kg]
ML_me = 4.48;  % Main engine mass and firing time [kg],[s]
t_me_L = 1763.9;

% Non recurring (development + one qualification unit)
NR_O_stm = 646*MO_stm^0.684;
NR_O_adcs = 324*MO_adcs;
NR_O_eps = 64.3*MO_eps;
NR_O_ps = 20*MO_rcs^0.485;
NR_O_tmtc = 618*MO_tmtc; % or 339*MO_tmtc + 5127*NO_ch;
NR_L_stm = 646*ML_stm^0.684;
NR_L_adcs = 324*ML_adcs;
NR_L_eps = 64.3*ML_eps;
NR_L_ps = 20*ML_rcs^0.485;
NR_L_tmtc = 618*ML_tmtc; % or 339*ML_tmtc + 5127*NL_ch;

NR_lander = (NR_L_stm+NR_L_adcs+NR_L_eps+NR_L_ps+NR_L_tmtc)*1000; 
NR_orbiter = (NR_O_stm+NR_O_adcs+NR_O_eps+NR_O_ps+NR_O_tmtc)*1000;

% Recurring (development of the flight unit)
R_O_stm = 22.6*MO_stm;
R_O_adcs = 795*MO_adcs^0.593;
R_O_eps = 32.4*MO_eps;
R_O_ps = 29*MO_me + t_me_O*0.024;
R_O_tmtc = 883.7*MO_tmtc^0.491;
R_L_stm = 22.6*ML_stm;
R_L_adcs = 795*ML_adcs^0.593;
R_L_eps = 32.4*ML_eps;
R_L_ps = 29*ML_me + t_me_L*0.024;
R_L_tmtc = 883.7*ML_tmtc^0.491;

R_lander = (R_L_stm+R_L_adcs+R_L_eps+R_L_tmtc)*1000; 
R_orbiter = (R_O_stm+R_O_adcs+R_O_eps+R_O_ps+R_O_tmtc)*1000;

% OBC unit cost
OBC_unit = 302000;   %[$], FY2023
OBC_O = 2*OBC_unit;
OBC_L = 3*OBC_unit;
OBC_cost = (OBC_O+OBC_L)*2; % doubled to take into account 1 QU and 1 FM

% FY conversion
NR_lander = NR_lander*1.2863; % $, FY2023
R_lander = R_lander*1.2863; % $, FY2023
NR_orbiter = NR_orbiter*1.2863; % $, FY2023
R_orbiter = R_orbiter*1.2863; % $, FY2023

R_ss = R_lander+R_orbiter;
NR_ss = NR_lander+NR_orbiter;
% Total hardware cost   $,FY2023
HW_cost = R_ss+NR_ss+PL_cost+OBC_cost;

%% Integration assembly and test [$, 2023]
NR_ait = 0.195*(NR_ss/1000+PL_cost/(2*1000)+OBC_cost/(2*1000)); 
R_ait = 0.124*(R_ss/1000+PL_cost/(2*1000)+OBC_cost/(2*1000)); 

C_iat = NR_ait*1000 + R_ait*1000; % [$, 2023]

%% Program level

NR_prog = 0.357*(NR_ss/1000+NR_ait/1000 +PL_cost/(2*1000));
R_prog = 0.320*(R_ss/1000+R_ait/1000 +PL_cost/(2*1000));

C_prog = NR_prog*1000+R_prog*1000; % [$, 2023]

%% Launcher cost

SLS = 2.2e+9;

%% Print results
fprintf('Orbiter Non-recurring cost: %d [M$,2023] \n',NR_orbiter/10^6);
fprintf('Lander Non-recurring cost: %d [M$,2023] \n',NR_lander/10^6);
fprintf('Orbiter recurring cost: %d [M$,2023] \n',R_orbiter/10^6);
fprintf('Lander recurring cost: %d [M$,2023] \n',R_lander/10^6);
fprintf('Total Non-recurring cost: %d [M$,2023] \n',NR_ss/10^6);
fprintf('Total recurring cost: %d [M$,2023] \n',R_ss/10^6);
fprintf('Total HW development cost: %d [M$,2023] \n',HW_cost/10^6);
fprintf('Total Program level cost: %d [M$,2023] \n',C_prog/10^6);
fprintf('Total IAT cost: %d [M$,2023] \n',C_iat/10^6);
fprintf('Main launcher service and vehicle cost: %d [M$,2023] \n',SLS/10^6);

%% Estimation error (SEE)
SEE_NR_stm = 0.22;
SEE_NR_adcs = 0.44;
SEE_NR_eps = 0.41;
SEE_NR_ps = 0.35;
SEE_NR_tmtc = 0.38;
SEE_NR_HW = sqrt(SEE_NR_stm^2+SEE_NR_eps^2+SEE_NR_ps^2+SEE_NR_adcs^2+...
    SEE_NR_tmtc^2);

SEE_R_stm = 0.21;
SEE_R_adcs = 0.36;
SEE_R_eps = 0.31;
SEE_R_ps = 0.22;
SEE_R_tmtc = 0.18;
SEE_R_HW = sqrt(SEE_R_stm^2+SEE_R_eps^2+SEE_R_ps^2+SEE_R_adcs^2+...
    SEE_R_tmtc^2);

SEE_NR_iat = 0.42; SEE_R_iat = 0.34;
SEE_iat = sqrt(SEE_NR_iat^2+SEE_R_iat^2);

SEE_NR_prog = 0.5; SEE_R_prog = 0.40;
SEE_prog = sqrt(SEE_NR_prog^2+SEE_R_prog^2);

