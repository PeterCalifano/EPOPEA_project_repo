clearvars; close all; clc

DefaultFontSize = 16;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

rng shuffle

%% Model parameters
J_SO = diag([12043.36, 15990.19, 8547.84]); % [kg m^2]

muS = 3.7931187*1e16;            %[m^3/s^2]
muE = (6.6743e-11)*(1.0802e20);  %[m^3/s^2]

% Geometry
inc = deg2rad(45);  % Thruster inclination wrt face 
lx = 3.600; % [m]
ly = 2.500; % [m]
lz = 3.610; % [m]
cm_off_x = 0.14; % [m]
cm_off_y = -0.01; % [m]
cm_off_z = -0.16; % [m]

D = -cm_off_z; % [m] thruster plane offset wrt CM CLAMPED (If not clamped put a plus)

J = J_SO;
ModelName = 'ADCS_orbiter';  

%% Main body normals
n1 = [1;0;0]; n4 = -n1; % x direction
n2 = [0;1;0]; n5 = -n2; % y direction
n3 = [0;0;1]; n6 = -n3; % z direction

% Main Body properties
A_cube = 12.9; % [m^2] considered taking 3.6 m x 3.6 m
rho_s = 0.5;
rho_d = 0.1;
r_cube = 0.5; % only valid for a 1x1x1 cube
% r_i = r_cube * [n1,n2,n3,n4,n5,n6];
r_i = [lx/2-cm_off_x  -cm_off_x       -cm_off_x        lx/2+cm_off_x  -cm_off_x       -cm_off_x
      -cm_off_y        ly/2-cm_off_y  -cm_off_y       -cm_off_y        ly/2+cm_off_y  -cm_off_y
      -cm_off_z       -cm_off_z        lz/2-cm_off_z  -cm_off_z       -cm_off_z        lz/2+cm_off_z ];

%% Initial conditions
wb0 = 0.1*rand(3, 1);
qrand = rand(4, 1);
qbn0 = qrand./norm(qrand);
% t0 = 0;
% tf = 1000;
% tspan = t0:0.1:tf;

% Initial interplanetary branch: from departure to first DSM
[R_PROP,V_PROP,TIMES,FLAGS, r_DSM, r_planets] = plot_trajectory(0);
% Adimensionalization variables
DU = astroConstants(2);
TU = 365.25;
TU2 = TU*3600*24;
VU = DU/TU2;
% Re-dimensionalize
R_PROP = R_PROP*DU;
V_PROP = V_PROP*VU;
r_DSM = r_DSM*DU;
r_planets = r_planets*DU;

% find when 3 h passed
hours = 3;
c = 0;
for j = 1:length(TIMES)
    check = TIMES(j) - TIMES(1);
    if check > (hours*3600)/TU2
        c = c + 1;
        final_time(c) = j; 
        break
    end
end

r_vec = R_PROP(:, 1:final_time);
v_vec = V_PROP(:, 1:final_time);
time = TIMES(1:final_time)*TU2;
TOF = time(end)-time(1);              % time of simulation
dt = time(2)-time(1);

% Interpolation to get vectors spanning 1s 
time_new = time(1):1:time(end);
dt_new = time_new(2)-time_new(1);
r_vec_new(1,:) = interp1(time,r_vec(1,:),time_new);
r_vec_new(2,:) = interp1(time,r_vec(2,:),time_new);
r_vec_new(3,:) = interp1(time,r_vec(3,:),time_new);
v_vec_new(1,:) = interp1(time,v_vec(1,:),time_new);
v_vec_new(2,:) = interp1(time,v_vec(2,:),time_new);
v_vec_new(3,:) = interp1(time,v_vec(3,:),time_new);

% from km to m and km/s to m/s
r_vec_new = r_vec_new*1e3;
v_vec_new = v_vec_new*1e3;

state0 = [r_vec_new; v_vec_new];

%% SRP Torque
n_sun = 2*pi/(365.25*24*3600);

eps = deg2rad(23.45);
ceps = cos(eps);
seps = sin(eps);

F_e = astroConstants(31); %W/m^2
c = 299792458; %m/s
P = F_e/c; %Pa

%% Sensors: Gyro only for detumbling
% Gyros model parameters
% Nsamples = length(tspan);
sigma_RRW = deg2rad(0.001);
sigma_ARW = deg2rad(0.001);
% R_RRW = diag(ones(1, 3)*sigma_RRW^2);
% R_ARW = diag(ones(1, 3)*sigma_ARW^2);
% RRW_noise = mvnrnd(zeros(3, 1), R_RRW, Nsamples);
% ARW_noise = mvnrnd(zeros(3, 1), R_ARW, Nsamples);

%% Actuators: Thrusters only for detumbling
% 8 actuators 

% MR 103J 1N Rocket Engine Thrusters
% MIB = 0.0133 Ns
F_min = 0.19;  % N
F_max = 1.13;  % N
MIB_real=F_min*0.1;
MIB_max=0.5;
Dt_min=MIB_max/F_min;

R = 0.5*[-ly*cos(inc)              ly*cos(inc)            ly*cos(inc)            -ly*cos(inc)             -ly*cos(inc)               ly*cos(inc)              ly*cos(inc)             -ly*cos(inc)          
          D*sin(inc)+lx*cos(inc)   D*sin(inc)-lx*cos(inc) D*sin(inc)+lx*cos(inc)  D*sin(inc)-lx*cos(inc)  -D*sin(inc)-lx*cos(inc)   -D*sin(inc)+lx*cos(inc)  -D*sin(inc)-lx*cos(inc)  -D*sin(inc)+lx*cos(inc)   
         -ly*sin(inc)             -ly*sin(inc)            ly*sin(inc)             ly*sin(inc)              ly*sin(inc)               ly*sin(inc)             -ly*sin(inc)             -ly*sin(inc)             ];

w = ones(length(R(1,:)),1);
R_ast = pinv(R);

%% Control 
% Detumbling requirement: 0.02 rad/s
% Gain scaling
G = 120;
k_dx = G*J(1,1)/max(diag(J));
k_dy = G*J(2,2)/max(diag(J));
k_dz = G*J(3,3)/max(diag(J));
k_p = 0;
k_d = 200;
w_c = [0 0 0]';

%% Run Simulink
fprintf('Model Ready! You can run Simulink')
% out = sim(ModelName, 'Timeout', StopTime);
