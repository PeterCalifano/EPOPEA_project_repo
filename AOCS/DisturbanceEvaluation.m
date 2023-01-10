close all;
clear;
clc;

% Masses: NSO, SL
Mass_arch = [7080.1; 7550.1]; 

J_NSO = diag([6597.09, 10138.44, 13254.21]); % [kg m^2]
J_SO = diag([12043.36, 15990.19, 8547.84]); % [kg m^2]
J = {J_NSO; J_SO};

cm_NSO = [-0.05; 1.51; -0.29]; %[m]
cm_SO = [-0.55; 1.75; 0.42]; %[m]
cm = [cm_NSO; cm_SO];

%% SRP

% Area
A_NS = [18.1, 16.705, 20.5];
A_S = [22.9, 21.875, 18];
A = [A_NS; A_S];

% Pressure center
disp_pc_cm = [3885.11, 4611.78]*1e-3;       % [m] NSO-SO

% Solar Radiation Pressure
AU = [1 9.6];                           % Earth - Saturn
c = 299792458;                          % [m/s]
I = 0;                                  % incidence angle
ref = 0.5;

T_SRPmax = zeros(2,2);
for j = 1:2               % rows location: Earth - Saturn
    for id = 1:2          % cols architecture: NSO - SO
        Area = max(A(id,:));
        SolarP = 1361*(1/AU(j)^2)/c;               % Pressure
        Fsrp = SolarP*Area*(1+ref)*cos(I);
        T_SRPmax(j,id) = Fsrp*disp_pc_cm(id);
    end
end

fprintf('\nSRP\n')
SRP = array2table([T_SRPmax(1,:); T_SRPmax(2,:) ],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Earth','Saturn'}); 
% Display table
disp(SRP) 

%% Magnetic Disturbance

% Magnetic Dipole
class = 1e-3;           % [Am^2/kg]
D = Mass_arch.*class;   % [Am^2]

% Magnetic field of Saturn at Enceladus distance (double check and report
% reference)
B_Enc = 325*1e-9;       %[T]
B_Earth = 40*1e-6;      %[T]
B = [B_Earth B_Enc];

for p = 1:2             % rows location: Earth - Saturn
    for i = 1:2         % cols architecture: NSO - SO
        T_Mmax(p,i) = D(i)*B(p);
    end
end

fprintf('\nMAGNETIC TORQUE\n')
M = array2table([T_Mmax(1,:); T_Mmax(2,:) ],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Interplanetary','Enceladus'}); 
disp(M)

%% Gravity Gradient

muS = 3.7931187*1e16;            %[m^3/s^2]
muE = (6.6743e-11)*(1.0802e20);  %[m^3/s^2]
muJ = 126686534	*1e9;            %[m^3/s^2]
mu = [muE muS muJ];

% Distance: minimum for Enceladus, mean from Saturn
min_d_E = 29000+251500;                 %[m]
mean_d_S = 238000000;            %[m] 
min_d_J = (69.911 + 2078300)*1e3;      %[m]
dist = [min_d_E mean_d_S min_d_J];

% Minimum and maximum inertia
Imin = [min(diag(J_NSO)) min(diag(J_SO))];
Imax = [max(diag(J_NSO)) max(diag(J_SO))];

% Tilting angle wrt Nadir: worst case
theta = deg2rad(45);

T_GGMax_div = zeros(2,2);
for p = 1:3         %rows: Enceladus, Saturn, Jupiter FlyBy
    for a = 1:2     %cols: NSO, SO
        T_GGMax_div(p,a) = 3*mu(p)/(2*dist(p)^3) * abs(Imin(a)-Imax(a)) * sin(2*theta);
    end
end
T_GGmax = [sum(T_GGMax_div(:,1)) sum(T_GGMax_div(:,2))];   % NSO - SO

fprintf('\nGRAVITY GRADIENT\n')
GG = array2table([T_GGMax_div(1,:); T_GGMax_div(2,:); T_GGMax_div(3,:); T_GGmax ],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Enceladus', 'Saturn', 'Jupiter flyby', 'Total @ Enceladus'});
disp(GG);

%% PLUME TORQUE
Cd = 2.1;
rho = 47*1e-12;             %[kg/m^3]
V = 213;                    %[m/s]
Across = [22.4 22.6];          %[m^2] NSO-SO
u_v = [0; 1; 0];
rr_cp_SO = [2.21; 2248.5; -1966]*1e-3;   %[m]
rr_cp_NSO = [-136; 1767.5; -1525]*1e-3;
rr_cp = [rr_cp_NSO rr_cp_SO];

for i=1:2
    TT_plume(i,:) = 0.5*Cd*rho*(V^2)*Across(i)*cross(u_v, rr_cp(:,i));
    T_plume(i) = norm(TT_plume(i,:));
end

% % Check if formula is correct
% T_cassini_true =  0.5*Cd*(5.49e-12)*(14410^2)*18.401*0.853;

fprintf('\nPLUME TORQUE\n')
total = array2table([T_plume],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Plume'});
disp(total);

%% TOTAL

% Interplanetary: SRP
Ttot_int = T_SRPmax(1,:) + T_Mmax(1,:);

% At Enceladus: GG + M +SRP
Ttot_Enc = T_SRPmax(2,:) + T_GGmax + T_Mmax(2,:) + T_plume;

fprintf('\nTOTAL TORQUE\n')
total = array2table([Ttot_int; Ttot_Enc],'VariableNames',{'NSO+SL','SO+SL'},'RowName',{'Interplanetary', '@ Enceladus'});
disp(total);

%% Actuator Sizing

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