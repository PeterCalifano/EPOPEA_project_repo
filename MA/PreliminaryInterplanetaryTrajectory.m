close all
clear
clc

%% Preliminary interplanetary Trajectory
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

% First guess: E-VEJ-S
% 1) Exit EARTH SOI with Vinf_FH (about 3 km/s) with transfer orbit Rp = 0.68 AU, Ra = 1 AU (slightly less)
% 2) GA VENUS with Vinf = about 5 km/s --> new transfer orbit Rp = 0.72 AU, Ra = about 1.02 AU
% 3) GA EARTH with Vinf = about 9 km/s --> new transfer orbit Rp = about 0.99 AU, Ra = about 5 AU

% 4) GA JUPITER with Vinf = about 6 km/s --> new transfer orbit Rp = about [3.6, 3.8] AU, Ra = about 9.9 AU
% 5) Entry SATURN SOI with Vinf = about 3 km/s --> Theoretical DV near 0 km/s if (3) is feasible
