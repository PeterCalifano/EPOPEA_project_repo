%%%%%%%%%%%%%%%%%%%%% STATION KEEPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load SPICE Kernels
cspice_kclear();
try
    cspice_furnsh('spice_kernels/pck00010.tpc')
    cspice_furnsh('spice_kernels/naif0012.tls')
    cspice_furnsh('spice_kernels/gm_de431.tpc')
    cspice_furnsh('spice_kernels/de440s.bsp')
    cspice_furnsh('spice_kernels/sat441.bsp')
catch
    cspice_furnsh('..\..\spice_kernels/pck00010.tpc')
    cspice_furnsh('..\..\spice_kernels/naif0012.tls')
    cspice_furnsh('..\..\spice_kernels/gm_de431.tpc')
    cspice_furnsh('..\..\spice_kernels/de440s.bsp')
end

%% DATA
clear; close all; clc;
% Constants
G = astroConstants(1);
mu_tbp = 1.90095713928102*1e-7;
DU = 238411468.296/1000; %km
TU = 118760.57/(2*pi); 

%mu_tbp = M_E/(M_E+M_S);

% Saturn and Enceladus Data
R_Saturn = astroConstants(26);
mu_Saturn = astroConstants(16);
R_Enceladus = mean(cspice_bodvrd('602','RADII',3));
mu_Enceladus = cspice_bodvrd('602','GM',1);
J2_Saturn = 1.629061510215236e-2; 
J2_Enceladus = 5435.2e-6; 
%mu_tbp = (mu_Enceladus)/(mu_Enceladus+mu_Saturn);

R_v = [R_Saturn, R_Enceladus]/DU;
mu_v = [mu_Saturn,mu_Enceladus] * TU^2 / DU^3;
J2_v = [J2_Saturn,J2_Enceladus];

%sample initial state for a resonant northern L2 orbit N=4, M=11
x0_Halo=1.000062853735440;
y0_Halo=0;
z0_Halo=-0.00117884381145460;
vx0_Halo=0;
vy0_Halo=0.0168877463349484;
vz0_Halo=0;
% x0_Halo = 1.000050323704250;
% z0_Halo = -0.00114508450748597;
% vy0_Halo = 0.0171721668270576;
% x0_Halo = 0.999937552716232;
% z0_Halo = -0.00118728429669647;
% vy0_Halo = -0.0168276813565369;

state0_Halo=[x0_Halo,y0_Halo,z0_Halo,vx0_Halo,vy0_Halo,vz0_Halo]';

t0=0;
FlightDays=0.5; %days of prapagation
tf=FlightDays*24*3600/TU; %final time of propagation
 
options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);

%propagation - Halo
[t_vec_Halo,state_vec_Halo]=ode113(@SCR3BP_dyn,[t0 tf],state0_Halo,...
    options_ode,mu_tbp,mu_v,R_v,J2_v);
state_vec_Halo=state_vec_Halo';

state_vec_Halo(1:3,:)=state_vec_Halo(1:3,:)*DU;
state_vec_Halo(4:6,:)=state_vec_Halo(4:6,:)*DU/TU;


Enceladus_3D(R_Enceladus,[(1-mu_tbp)*DU,0,0]);
P2=plot3(state_vec_Halo(1,:),state_vec_Halo(2,:),state_vec_Halo(3,:),...
    'k','linewidth',1.25,'DisplayName','Southern Halo Orbit');
%S = scatter3(-mu_tbp*DU,0,0,100,'filled','DisplayName','Saturn');
%P1=plot3(x_L2(1)*DU,x_L2(2)*DU,x_L2(3)*DU,'ob','markersize',5,'linewidth',1.25);
grid minor

%% STATION KEEPING
N = 4 % Select the number of control points for the SK during one orbit
states_SK = zeros(6,N); % Adimensional SK states
%times_SK = zeros(1,N); % Adimensional SK times


%%% INFORMATION ABOUT THE DIFFERENT ARCS
h_RS = 1; %[h] - duration of the remote sensing arc (1 of the 3 modes)
tf_RS=h_RS/2*3600/TU; 
h_CI=2; %[h] - duration of the coarse imaging arc
tf_CI=tf_RS+h_CI/2*3600/TU; 
h_SK=2; %Number of hours dedicated to SK
tf_SK = tf_CI+h_SK/2*3600/TU;
t_half = 1.142397328535602;
ti_SK2 = t_half + (t_half - tf_SK);
tf_SK2 = 2*t_half - tf_CI;
times_SK = [tf_CI,tf_SK,ti_SK2,tf_SK2];

% Define integration options
options_ode_event = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13,'Events',@(t,x) FirstZeroCrossing(t,x));
options_ode = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13);

% Find nominal states in the SK Points
t0_prop_ii = t0;
state0_prop_ii = state0_Halo;
for ii = 1:N

    % Propagate until time of interest
    [t_ii,state_v_ii] = ode113(@(t,x) CR3BP_dyn(t,x,mu_tbp),[t0_prop_ii, times_SK(ii)],state0_prop_ii, options_ode);

    % Save the state
    states_SK(:,ii) = state_v_ii(end,:)';

    % Update for the next propagation
    t0_prop_ii = times_SK(ii);
    state0_prop_ii = states_SK(:,ii);
end

% PERFORM SK

for ii = 1:N
    



end


% Propagation of a full nominal Halo
[t_v,state_fullHalo]=ode113(@(t,x) CR3BP_dyn(t,x,mu_tbp),[t0, 2*t_half],state0_Halo, options_ode);
state_fullHalo(:,1:3) = state_fullHalo(:,1:3)*DU;
state_fullHalo(:,4:6) = state_fullHalo(:,4:6)*DU/TU;

% Plot
Enceladus_3D(R_Enceladus,[(1-mu_tbp)*DU,0,0]);
P2=plot3(state_fullHalo(:,1),state_fullHalo(:,2),state_fullHalo(:,3),...
    'k','linewidth',1.25,'DisplayName','Nominal Halo Orbit');
grid minor
hold on
scatter3(states_SK(1,:)*DU,states_SK(2,:)*DU,states_SK(3,:)*DU,40,'r','filled',...
    'DisplayName','SK Points')
legend()



function [value,isterminal,direction] = FirstZeroCrossing(t,x) %Event function 
%to stop the integration at the intersection point with the x-z plane 

value=x(2);
isterminal=1;
direction=-1;

end
