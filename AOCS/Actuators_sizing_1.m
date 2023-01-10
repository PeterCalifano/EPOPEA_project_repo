%% Actuator Sizing

% Mass
Mass_arch = 7550.1; %%%%%%%%%%%%% CHANGE %%%%%%%%%%%%%%
% Inertia
J_SO = diag([6846.3, 15903.78, 18425.41]); % [kg m^2]
% Center of mass position with respect to the base
cm_SO = [-0.55; 1.75; 0.42]; %[m]

%% Reaction Wheel

% Slew
% RATE: 0.05 dps
rotation = deg2rad(30);            %[rad]
time = 600;                        %[s]

for i = 1:2          % NS - S
    Ji = max(diag(J{i,:}));
    T_RW_slw(i) = 4*rotation*Ji/time^2;
end

% Momentum storage
orb_period = 12*3600;

for i = 1:2
    H(i) = T_GGmax(i)*orb_period/4;
end

%% Thrusters

% Slew
rot = deg2rad(30);                 %[rad]
time = 60;                         %[s]
rate = rot/time;                   %[rad/s]
time_acc = 0.05*time;              %[s] 
acc = rate/time_acc;

b = [3500 2960          %NSO
     3600 2500].*1e-3;        %SO
for i = 1:2
    Ji = max(diag(J{i,:}));
    T = Ji*acc;

    % Two thrusters with highest momentum arm
    F_2T(i) = T/b(i,1);

    % Four thrusters with F1 = 0.75 F2 (L1 > L2)
    F_2(i) = T/(b(i,2)+0.75*b(i,1));
    F_1(i) = 0.75*F_2(i);
end