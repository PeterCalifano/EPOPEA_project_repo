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
DU=238411468.296/1000; %km
TU=118760.57/(2*pi); 

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
FlightDays=5; %days of prapagation
tf=FlightDays*24*3600/TU; %final time of propagation
 
options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);

%propagation - Halo
[t_vec_Halo,state_vec_Halo]=ode113(@CR3BP_dyn,[t0 tf],state0_Halo,...
    options_ode,mu_tbp);%,mu_v,R_v,J2_v);
state_vec_Halo=state_vec_Halo';

state_vec_Halo(1:3,:)=state_vec_Halo(1:3,:)*DU;
state_vec_Halo(4:6,:)=state_vec_Halo(4:6,:)*DU/TU;


Enceladus_3D(R_Enceladus,[(1-mu_tbp)*DU,0,0])
P2=plot3(state_vec_Halo(1,:),state_vec_Halo(2,:),state_vec_Halo(3,:),...
    'k','linewidth',1.25,'DisplayName','Southern Halo Orbit');
%S = scatter3(-mu_tbp*DU,0,0,100,'filled','DisplayName','Saturn');
%P1=plot3(x_L2(1)*DU,x_L2(2)*DU,x_L2(3)*DU,'ob','markersize',5,'linewidth',1.25);
grid minor

%% STATION KEEPING
N = 4 % Select the number of control points for the SK during one orbit
states_SK0 = zeros(6,N); % Adimensional SK states
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
times_SK0 = [tf_CI,tf_SK,ti_SK2,tf_SK2];

% Define integration options
options_ode_event = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13,'Events',@(t,x) FirstZeroCrossing(t,x));
options_ode = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13);

% Find nominal states in the SK Points
t0_prop_ii = t0;
state0_prop_ii = state0_Halo;
for ii = 1:N

    % Propagate until time of interest
    [t_ii,state_v_ii] = ode113(@(t,x) CR3BP_dyn(t,x,mu_tbp),[t0_prop_ii, times_SK0(ii)],state0_prop_ii, options_ode);

    % Save the state
    states_SK0(:,ii) = state_v_ii(end,:)';

    % Update for the next propagation
    t0_prop_ii = times_SK0(ii);
    state0_prop_ii = states_SK0(:,ii);
end

%% SIMPLE SHOOTING
clc
options = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'iter',...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6, 'ConstraintTolerance', 1e-6,...
    'SpecifyObjectiveGradient', false, 'SpecifyConstraintGradient', false, ...
    'MaxFunctionEvaluations',5000000,'MaxIterations',500000,'FunctionTolerance',1e-6); 

N_orbits = 1;

states_SK = [];
times_SK = [];

for i = 1 : N_orbits
    states_SK = [states_SK, states_SK0];
    times_SK = [times_SK, i * 2*t_half + times_SK0];
end

n_var = 8; % 6 components of state + 2 times

N = length(times_SK);

% Set thresholds
threshold_t = 5 * 60/TU;
threshold_r = 100/DU; 

A = [];
B = [];
Aeq = [];
Beq = [];
ub = 1e+10*ones(n_var,1);
lb = -1e+10*ones(n_var,1);

state_output = zeros(n_var-2,N-1);
times_output = zeros(2,N-1);
DV_output = zeros(1,N-1);

state0 = states_SK(:,1);
time0 = times_SK(1);

for ii = 1: N-1

    if rem(ii,4) == 0
        flag = 1; % Pericenter
        bounds = [20,60]/DU;
    elseif rem(ii,2) == 0 && rem(ii,4) ~= 0
        flag = 2; % Apocenter
        bounds = [800,1500]/DU;
    else
        flag = 0; % SK arcs
        bounds = [0,0];
    end

    r_nom_i = states_SK(1:3,ii);
    r_nom_f = states_SK(1:3,ii+1);
    v_nom_i = states_SK(4:6,ii);

    x0 = [state0 ;
        time0;
        times_SK(ii+1) - times_SK(ii) + time0];
    lb(7) = times_SK(ii) - threshold_t;
    ub(7) = times_SK(ii) + threshold_t;
    lb(8) = times_SK(ii + 1) - threshold_t;
    ub(8) = times_SK(ii + 1) + threshold_t;

    [X_ss, DV_ss] = fmincon(@(var) objfun_SK(var,mu_tbp,mu_v,R_v,J2_v,v_nom_i),x0,A,B,...
        Aeq,Beq,lb,ub,@(var) nlcon_SK(var,mu_tbp,mu_v,R_v,J2_v,r_nom_i,r_nom_f,threshold_r,flag,bounds),options);
    
    state_output(:,ii) = X_ss(1:6);
    times_output(1,ii) = X_ss(7);
    times_output(2,ii) = X_ss(8);
    
    DV_output(ii) = DV_ss;
    
    state0 = state_output(:,ii);
    time0 = times_output(2,ii);
end

DV_output = DV_output * DU/TU * 1000;


%% Propagation of a full nominal Halo
[t_v,state_fullHalo]=ode113(@SCR3BP_dyn,[t0, 2*t_half],state0_Halo, options_ode...
    ,mu_tbp,mu_v,R_v,J2_v);
state_fullHalo(:,1:3) = state_fullHalo(:,1:3)*DU;
state_fullHalo(:,4:6) = state_fullHalo(:,4:6)*DU/TU;


% Plot
Enceladus_3D(R_Enceladus,[(1-mu_tbp)*DU,0,0]);
P2=plot3(state_fullHalo(:,1),state_fullHalo(:,2),state_fullHalo(:,3),...
    'k--','linewidth',0.5,'DisplayName','Nominal Halo Orbit');
grid minor
hold on
scatter3(states_SK0(1,:)*DU,states_SK0(2,:)*DU,states_SK0(3,:)*DU,40,'r','filled',...
    'DisplayName','SK Points')
scatter3(state_output(1,:)*DU,state_output(2,:)*DU,state_output(3,:)*DU,40,'g','filled',...
    'DisplayName','Optim 1')


for ii = 1:N-1
    xx_0 = state_output(:,ii);
    t0 = times_output(1,ii);
    tf = times_output(2,ii);
    [~,prop_state] = ode113(@SCR3BP_dyn,[t0, tf],xx_0,options_ode,mu_tbp,mu_v,R_v,J2_v);
    prop_state = prop_state';
    plot3(prop_state(1,:)*DU,prop_state(2,:)*DU,prop_state(3,:)*DU,'DisplayName','Propagation')
end
legend()

%% MULTIPLE SHOOTING
clc
options = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'iter',...
    'OptimalityTolerance', 1e-9, 'StepTolerance', 1e-9, 'ConstraintTolerance', 1e-9,...
    'SpecifyObjectiveGradient', false, 'SpecifyConstraintGradient', false, ...
    'MaxFunctionEvaluations',5000000,'MaxIterations',500000,'FunctionTolerance',1e-9); 

% var [14*N_orbits]

lb_peri = 20/DU;
ub_peri = 60/DU;
lb_apo = 800/DU;
ub_apo = 1500/DU;

N_orbits = 1

n_var = 14*N_orbits;

t_orb = 2*t_half;

% Create linear constraints and boundaries
A = [];
B = [];
Aeq = [];
Beq = [];
ub = 1e+10*ones(n_var,1);
lb = -1e+10*ones(n_var,1);

for ii = 1:2*N_orbits
    
    if rem(ii,2) ~= 0

        % SK1 ( peri --> apo )
        ub(7*ii) = (ii-1)/2 * t_orb + tf_SK;
        lb(7*ii) = (ii-1)/2 * t_orb + tf_CI;
        lb(7*ii - 5) = 0;

    else

        % SK2 ( apo --> peri )
        ub(7*ii) = ii/2 * t_orb - tf_CI;
        lb(7*ii) = ii/2 * t_orb - tf_SK;
        ub(7*ii - 5) = 0;

    end

end

% Create bounds for SK position and velocity
norm_r1 = norm(states_SK0(1:3,1) - [1-mu_tbp;0;0]);
norm_r2 = norm(states_SK0(1:3,2) - [1-mu_tbp;0;0]);
norm_v1 = norm(states_SK0(4:6,1));
norm_v2 = norm(states_SK0(4:6,2));

r_max = max(norm_r1,norm_r2)*1.2;
r_min = min(norm_r1,norm_r2)*0.8;
v_max = max(norm_v1,norm_v2)*1.2;
v_min = min(norm_v1,norm_v2)*0.8;

% Create initial guess
initial_guess = zeros(n_var,1);

for ii = 1:2*N_orbits

    if rem(ii,2) ~= 0
        index = 2;
        sum = -0.001;
    else
        index = 3;
        sum = 0.001;
    end

    initial_guess(7*ii-6:7*ii-1) = states_SK0(:,index);
    initial_guess(7*ii) = times_SK0(index) + sum;

end

[X_ms, DV_ms] = fmincon(@(var) objfcn_multiple_SK(var,mu_tbp,mu_v,R_v,J2_v,state0_Halo,N_orbits),initial_guess,A,B,...
    Aeq,Beq,lb,ub,@(var) nlcon_multiple_SK(var,mu_tbp,mu_v,R_v,J2_v,state0_Halo,N_orbits,lb_peri,ub_peri,lb_apo,ub_apo,...
    r_max,r_min,v_max,v_min), options);
    



