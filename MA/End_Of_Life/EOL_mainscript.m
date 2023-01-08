%%%%%%%%%%%%%%%%%%%%% END OF LIFE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);
options_ode_event = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13,'Events',@FirstZeroCrossing);

% General data
R_Enceladus = mean(cspice_bodvrd('602','RADII',3));
R_Saturn = mean(cspice_bodvrd('699','RADII',3));
mu_Saturn = cspice_bodvrd('699','GM',1);
mu_Enceladus = cspice_bodvrd('602','GM',1);
mu_tbp = 1.90095713928102*1e-7;
J2_Saturn = 1.629061510215236e-2; 
J2_Enceladus = 5435.2e-6; 
DU=238411468.296/1000;
TU=118760.57/(2*pi);  
VU = DU/TU;
t_half = 1.142397328535602;
R_v = [R_Saturn, R_Enceladus]/DU;
mu_v = [mu_Saturn,mu_Enceladus] * TU^2 / DU^3;
J2_v = [J2_Saturn,J2_Enceladus];

% Initial Halo State
x0_Halo=1.000062853735440;
y0_Halo=0;
z0_Halo=-0.00117884381145460;
vx0_Halo=0;
vy0_Halo=0.0168877463349484;
vz0_Halo=0;
state0_Halo=[x0_Halo,y0_Halo,z0_Halo,vx0_Halo,vy0_Halo,vz0_Halo]';
t0 = 0;

% Design variables
% var = [x_1, t_1, t_f];
n_var = 15;

% Set the bounds
A = zeros(3,n_var);
A(1,[7,14]) = [1,-1];
A(2,[14,15]) = [1, -1];
A(3,7) = -1;
B = [-0.1,-0.1,-t0-0.1];
Aeq = [];
Beq = [];
ub = 1e+10*ones(n_var,1);
lb = -1e+10*ones(n_var,1);
lb(7) = 0;


% Create Initial Guess
[~,prop_state] = ode113(@CR3BP_dyn,[t0 t_half-0.2],state0_Halo,options_ode,mu_tbp);
x1_guess = prop_state(end,:)';

max_biell = 2.6;
x2_in = [-max_biell*R_Saturn,0,0,0,-sqrt(mu_Saturn/(max_biell*R_Saturn)),0];
a = (max_biell*R_Saturn + 4*R_Saturn)/2;
t2 = 2*pi*sqrt(a^3/mu_Saturn) / 2 / TU;
x2_guess(1) = + (x2_in(1)-mu_tbp)*cos(t2) + x2_in(2)*sin(t2);
x2_guess(2) = - (x2_in(1)-mu_tbp)*sin(t2) + x2_in(2)*cos(t2);
x2_guess(3) = x2_in(3);
x2_guess(4) = x2_in(4);
x2_guess(5) = x2_in(5);
x2_guess(6) = x2_in(6);

x2_guess(1:3) = x2_guess(1:3)/DU;
x2_guess(4:6) = x2_guess(4:6)*TU/DU;

%x2_guess = [1.00195813111270, 0.000330552288712704 0.00484092359130717, ...
%    -0.00124257763903574, -0.110815788777326, -0.000139707519341356];
%x2_guess = [0.999961968882375,-5.61855695907048e-05, -6.22774331811784e-05...
%    -0.00980074719104070, -0.122875589208661, -0.00140440124794177];

a = (R_Saturn + 4*R_Saturn)/2;
T = 2*pi*sqrt(a^3/mu_Saturn);
tf_guess = T / TU;

initial_guess = [x1_guess; t_half - 0.3; x2_guess'; 2*t_half; 4 * t_half];

load('InitialGuess_EOL_DV1.34.mat','XX_eol')
initial_guess = XX_eol;
R_final = 2.7*R_Saturn/DU;
r_max = 300*R_Saturn/DU;

%% Optimization
options = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'iter',...
    'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, 'ConstraintTolerance', 1e-11,...
    'SpecifyObjectiveGradient', false, 'SpecifyConstraintGradient', false, ...
    'MaxFunctionEvaluations',5000000,'MaxIterations',500000,'FunctionTolerance',1e-8); 

[XX_eol, DV_eol] = fmincon(@(var) objfun_EOL(var,mu_tbp,t0,state0_Halo,mu_v,R_v,J2_v,DU/TU),initial_guess,A,B,...
    Aeq,Beq,lb,ub,@(var) nlcon_EOL(var,mu_tbp,t0,state0_Halo,R_final,r_max,mu_v,R_v,J2_v), options);


DV_opt = objfun_EOL(XX_eol,mu_tbp,t0,state0_Halo,mu_v,R_v,J2_v,DU/TU)

%% Post-Processing
x1_opt = XX_eol(1:6);
t1_opt = XX_eol(7);
x2_opt = XX_eol(8:13);
t2_opt = XX_eol(14);
tf_opt = XX_eol(15);

[t_prop0,zero_arc] = ode113(@CR3BP_dyn,[t0, t1_opt],state0_Halo,options_ode,mu_tbp);
[t_prop1,first_arc] = ode113(@SCR3BP_dyn,[t1_opt t2_opt],x1_opt,options_ode,mu_tbp,mu_v,R_v,J2_v);
[t_prop2,second_arc] = ode113(@SCR3BP_dyn,[t2_opt tf_opt],x2_opt,options_ode,mu_tbp,mu_v,R_v,J2_v);
[t_propf,final_arc] = ode113(@SCR3BP_dyn,[tf_opt, 2*tf_opt],second_arc(end,:)',options_ode,mu_tbp,mu_v,R_v,J2_v);
prop_state = [zero_arc',first_arc',second_arc',final_arc'];
t_prop = [t_prop0; t_prop1;t_prop2;t_propf];

fa = length(zero_arc(:,1)) + length(first_arc(:,1));
fa2 = length(zero_arc(:,1)) + length(first_arc(:,1)) + length(first_arc(:,2));

prop_state_in = zeros(3,length(t_prop));
for k=1:length(t_prop)
    prop_state_in(1,k)=(prop_state(1,k)+mu_tbp)*cos(t_prop(k) -t_prop(1))-prop_state(2,k)*sin(t_prop(k) -t_prop(1));
    prop_state_in(2,k)=(prop_state(1,k)+mu_tbp)*sin(t_prop(k) -t_prop(1))+prop_state(2,k)*cos(t_prop(k) -t_prop(1));
    prop_state_in(3,k)=prop_state(3,k);
end

Enceladus_3D(R_Saturn,[0,0,0]);
P1 = plot3(prop_state_in(1,1:fa)*DU,prop_state_in(2,1:fa)*DU,prop_state_in(3,1:fa)*DU,...
    'b','linewidth',2);
P2 = plot3(prop_state_in(1,fa:fa2)*DU,prop_state_in(2,fa:fa2)*DU,prop_state_in(3,fa:fa2)*DU,...
    'r','linewidth',2);
P3 = plot3(prop_state_in(1,fa2:end)*DU,prop_state_in(2,fa2:end)*DU,prop_state_in(3,fa2:end)*DU,...
    'g--','linewidth',2);
%scatter3((1-mu_tbp)*DU,0,0,50,'k','filled','DisplayName','Enceladus')

ind0 = length(zero_arc(:,1));
ind1 = length(zero_arc(:,1)) + 50;

Enceladus_3D(R_Enceladus,[(1-mu_tbp)*DU,0,0]);
P3 = plot3(prop_state(1,1:ind0)*DU,prop_state(2,1:ind0)*DU,prop_state(3,1:ind0)*DU,...
    'b','linewidth',2);
P4 = plot3(prop_state(1,ind0:ind1)*DU,prop_state(2,ind0:ind1)*DU,prop_state(3,ind0:ind1)*DU,...
    'r','linewidth',2);

%% 
N_targ = 3.5;

A2 = [];
B2 = [];
lb2 = 2;
ub2 = 200;
[N_opt, DV_opt] = fmincon(@(N) biell_cost(N,N_targ,mu_v,R_v,DU/TU),20, A2,B2,...
    Aeq,Beq,lb2,ub2,[], options);

DV_opt


