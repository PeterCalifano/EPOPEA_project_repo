%% Actuator Sizing

% Mass
% Mass = 7550.1; %%%%%%%%%%%%% CHANGE %%%%%%%%%%%%%%
% Inertia Clamped 
J = diag([12043.36, 15990.19, 8547.84]); % [kg m^2]
% Center of mass position with respect to the base (adapter face)
cm_SO = [0.14; -0.01; 1.64]; %[m]


%% Thrusters
% Formula: I * 𝜃.. = Fd

inc = deg2rad(45);  % Thruster inclination wrt face 
lx = 3.600; % [m]
ly = 2.500; % [m]
lz = 3.610; % [m]

F_nom = 1; % [N] nominal thrust range: 0.19 N -1-13 N
acc_time = 5; % [s] acceleration pulse
rot_y = deg2rad(90);
rot_x = deg2rad(90);
rot_z = deg2rad(180);

% rotation around y
d_y = lx*cos(inc);
dd_theta_y = (4*F_nom*d_y)/J(2,2);
d_theta_y = dd_theta_y * acc_time;
time_rot_y = rot_y/d_theta_y/60; % [min] time to perform a rotation equal to rot_y

% rotation around x
d_x = ly;
dd_theta_x = (4*F_nom*cos(inc)*d_x)/J(1,1);
d_theta_x = dd_theta_x * acc_time;
time_rot_x = rot_x/d_theta_x/60; % [min] time to perform a rotation equal to rot_y

% rotation around z
d_z = ly;
dd_theta_z = (4*F_nom*sin(inc)*d_z)/J(3,3);
d_theta_z = dd_theta_z * acc_time;
time_rot_z = rot_z/d_theta_z/60; % [min] time to perform a rotation equal to rot_y

%% Thruster pulse life

