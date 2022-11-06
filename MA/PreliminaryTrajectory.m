close all
clear
clc


cspice_kclear();
cspice_furnsh('..\..\EPOPEA_project_repo\EPOPEA_metakernel.tm')


%% Preliminary interplanetary Trajectory
mu_Sun = cspice_bodvrd('SUN', 'GM', 1); % [km/s^2]

m0_assumed = 5000; % [kg]
C3_FH = 11; % [km^2/s^2]
Vinf_FH = sqrt(C3_FH); % [km/s]

% GRAPHS Rp-Ra ARE CREATED WITH ALPHA = 180° at the bottom end of the curves (0° at the top)
% ASSUMPTION: PLANAR AND CIRCULAR ORBIT
% First guess: E-VEJ-S
% 1) Exit EARTH SOI with Vinf_FH (about 3 km/s) with transfer orbit Rp = 0.68 AU, Ra = 1 AU (slightly less)
% 2) GA VENUS with Vinf = about 5 km/s --> new transfer orbit Rp = 0.72 AU, Ra = about 1.02 AU
% 3) GA EARTH with Vinf = about 9 km/s --> new transfer orbit Rp = about 0.99 AU, Ra = about 5 AU --> feasibility to be assessed (max alpha constraint)
% 4) GA JUPITER with Vinf = about 6 km/s --> new transfer orbit Rp = about [3.6, 3.8] AU, Ra = about 9.9 AU
% 5) Entry SATURN SOI with Vinf = about 3 km/s --> Theoretical DV near 0 km/s if (3) is feasible

% Second guess: E-VEVE-S
% 1) Exit EARTH SOI with Vinf_FH (about 3 km/s) with transfer orbit Rp =  0.68 AU, Ra = 1 AU 
% 2) GA VENUS Vinf = 5 km/s --> Ra = 1.3 AU, Rp = 0.71 AU
% 3) GA EARTH Vinf = 8 km/s --> Ra = 1.2 AU, Rp = 0.67 AU
% 4) GA VENUS Vinf = 7 km/s --> Ra = 1.5 AU, Rp = 0.7 AU
% 5) GA EARTH Vinf = 11 km/s --> Ra = 9.1 AU, Rp = 0.98 AU
% 6) Entry SATURN SOI with Vinf = about 6 km/s

% Third guess: E-VEE-S
% 1) Exit EARTH SOI with Vinf_FH (about 4 km/s) with transfer orbit Rp = 0.61 AU, Ra = 1 AU 
% 2) GA VENUS Vinf = 8 km/s --> Ra =  1.5 AU, Rp = 0.7 AU
% 3) GA EARTH Vinf = 11 km/s --> Ra =  2.4 AU, Rp = 0.86 AU
% 5) GA EARTH Vinf = 11 km/s --> Ra = 9.1 AU, Rp = 0.98 AU
% 6) Entry SATURN SOI with Vinf = about 6 km/s


%% Preliminary Saturn System trajectory
mu_Saturn = cspice_bodvrd('SATURN', 'GM', 1); % [km/s^2]

R_Saturn = mean(cspice_bodvrd('SATURN', 'RADII', 3)); % [km]
R_Enceladus = mean(cspice_bodvrd('602', 'RADII', 3)); % [km]

Ra_target = 200*R_Saturn;
Rp_target = 20*R_Saturn;
Vinf_arrival = 3; % [km/s]

% Hypothesis: energy evaluation only
% Parameters of target orbit
SMA_target = (Ra_target + Rp_target)/2;
eps_target = -mu_Saturn/(2*SMA_target);
e_target = (Ra_target - Rp_target)/(Ra_target + Rp_target);

% Velocity of hyperbolic orbit at Rp target
Vp_hyp = sqrt(Vinf_arrival^2 + 2*mu_Saturn/(Rp_target)); 

% Velocity at pericentre of target orbit
V_cap = sqrt(mu_Saturn*(1 + e_target)/Rp_target);

% DV required to achieve elliptical capture orbit
dV_capture = Vp_hyp - V_cap

% Period of the capture orbit
T_target = 2*pi*sqrt(SMA_target^3./mu_Saturn)/3600; % [h]
T_target_d = T_target./24; % [days]

% transfer_method = 1;
% transfer_arc_switch
% [deltaV_arrival, deltaV_injection, deltaV_tot, transfer_arc] = orbit_transfer(mu, r_i1, r_f2, v_i1, v_f2, ToF, transfer_method, transfer_arc_switch);







