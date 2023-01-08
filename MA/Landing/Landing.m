%% Landing problem: Collocation method
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(0,'defaultAxesFontSize', 16)

%%
clearvars; close all; 
clc 

n_T = 4; % number of thrusters.  try with 3 thrusters
% step1: vehicle par (TO CHANGE!!!)
Tmax = n_T*220e-3; %[kN] Maximum Thrust  !!!!! 4* 
Isp = 229; %[s] Specific Impulse
g0 = 9.81*1e-3; %[km/s^2] acceleration constant

% NEW MASS:
m_dry_lander =508; % [kg]
m_prop_HA = 80; % [kg] mass of propellant for hazard avoidance. NOTE: HYP
m_dry = m_dry_lander + m_prop_HA; % NOTE: improper use of name "dry mass": 
% it is the lower bound for the mass of the lander. 
% We cannot use the propellant dedicated to hazard avoidance maneouvre to land
m_prop = 252; % [kg] NOTE: estimated mass of propellant from PS, can be iterated
m0 = m_dry_lander+ m_prop; % [kg] WET MASS: initial mass of landing trajectory

% Enceladus par
Re = 251.1; %[km] mean radius of Enceladus
mu = (6.67430e-11 *1.0802e20 )*10^(-9) ; % [km^3/s^2] Enceladus gravitational constant
mass_ratio = 1.90095713928102*1e-7; % Saturn-Enceladus mass ratio
we = (1/32.9)*2*pi/3600; %[rad/s] enceladus angular rate

% Adimensionalization
DU = Re;                        %[km]
TU = 3600;                      %[s]
VU = DU/TU;                     %[km/s]
MM = m0;                        %[kg]
acc = DU/(TU)^2;                %[km/s^2]
FU = MM*acc;                    %[kN]
GM = DU^3/TU^2;                 %[km^3/s^2]

Tmax = Tmax/FU;
Isp = Isp/TU;
g0 = g0/acc;
m0 = m0/MM;
m_dry = m_dry/MM;
Re = Re/DU;
mu = mu/GM;
we = we*TU;

% organize par 
par(1) = Tmax;
par(2) = Isp;
par(3) = g0;
par(4) = mu;
par(5) = Re;
par(6) = mass_ratio;
par(7) = we;

N = 20;

% INITIAL ORBIT : circular, polar (TO CHANGE!!)
initial = 3;
switch initial
    case 1
        h = 60/DU;
        r_mod = Re+h;
        xi = -r_mod/sqrt(2);
        yi = -r_mod/sqrt(2);
        zi = 0;
        v_mod = sqrt(mu/norm([xi; yi]));
        vxi = v_mod/sqrt(2);
        vyi = -v_mod/sqrt(2);
        vzi = 0;
        state_i = [xi yi zi vxi vyi vzi m0];
        tN = 0.7*3600/TU;              %final time guess if h = 100 km
        % Initial time guess
        t1 = eps;            
    case 2
        DD=238411468.296/1000;       %km
        TT=118760.57/(2*pi);         %s
        xi = 238479.035029632/DD;
        yi = -443.388268197043/DD;
        zi = -57.8313050277481/DD;
        vxi = -0.0330981241105157/(DD/TT);
        vyi = 0.103657810950128/(DD/TT);
        vzi = -0.120139871111334/(DD/TT);
        state_i_rot = [xi yi zi vxi vyi vzi m0];
        % Initial time guess
        t1 = eps; 
        % Final time guess 
        tN = 0.7 + t1;   

        % From rotating Saturn-Enceladus to IAU_Enceladus
        state_in = rot2iau_enc(t1*(TU/TT), state_i_rot(1:end-1), mass_ratio);
        r_in = state_in(1:3)*DD/DU;
        v_in = state_in(4:6)*(DD/TT)/VU;
        
        % DCM Matrix: rigid rotation around X
        A_rotx = [1  0  0
                  0  0  1
                  0 -1  0];
        r1 = A_rotx*r_in;
        v1 = A_rotx*v_in;
        state_i = [r1' v1' m0];
    case 3
        DD=238411468.296/1000;       %km
        TT=118760.57/(2*pi);         %s
        xi = 238499.809917516/DD;
        yi = -4.989057360755955e+02/DD;
        zi = 14.352712909418244/DD;
        vxi = -0.036329671799006/(DD/TT);
        vyi = 0.082968629856578/(DD/TT);
        vzi = -0.121113200205420/(DD/TT);
        state_i_rot = [xi yi zi vxi vyi vzi m0];
        % Initial time guess
        t1 = eps; 
        % Final time guess 
        tN = 0.7 + t1;   

        % From rotating Saturn-Enceladus to IAU_Enceladus
        state_in = rot2iau_enc(t1*(TU/TT), state_i_rot(1:end-1), mass_ratio);
        r_in = state_in(1:3)*DD/DU;
        v_in = state_in(4:6)*(DD/TT)/VU;
        
        % DCM Matrix: rigid rotation around X
        A_rotx = [1  0  0
                  0  0  1
                  0 -1  0];
        r1 = A_rotx*r_in;
        v1 = A_rotx*v_in;
        state_i = [r1' v1' m0];
end

% Landing site
landing_site = 3;
switch landing_sit
    case 1
        lonlat = [-80; 20]; 
        th_e0 = 0;
    case 2
        lonlat = [-70; 270]; 
        th_e0 = deg2rad(290);
    case 3
        lonlat = [-70; 270]; 
        th_e0 =deg2rad(270);
    case 4
        t1 = -0.2;
        lonlat = [-68;20];  
        th_e0 =deg2rad(-10);
    case 5
        lonlat = [-90; 0]; 
        th_e0 = 0;
end
par(8) = th_e0;

% NLP vars (x1, u1, ..., xN, uN, t1, tN)
step_st = length(state_i);            % 7: rr,vv,m
step_var = step_st + 4;               % 11: rr,vv,m,u,ax,ay,az
                            

% Time grid and initial guess fill
h = (tN-t1)/(N-1);
s0 = state_i;
guess = zeros(step_var*N+2, 1);
guess(1:step_st,1) = s0;
guess(step_st+1) = 1-eps;
guess(step_st+2:step_st+4) = -s0(4:6)'/norm(s0(4:6));
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

tspan_l(1) = t1;
tspan_l(N) = tN;
for k = 2:N-1
    tspan_l(k) = t1 + h*(k-1);
end

uk = guess(step_st+1:step_st+4);
for k = 1:N-1
    [~, output] = ode113(@landing_dyn, [tspan_l(k) tspan_l(k+1)], s0', options, uk, par);
    s0 = output(end,1:step_st);
    vv_vers = s0(4:6)'/norm(s0(4:6));
    if k == N-1
        uk = [1-eps; -vv_vers];
    else
        uk = [eps; -vv_vers];
    end
    guess(k*step_var+1:(k+1)*step_var) = [s0'; uk];
end
guess(end-1:end) = [t1; tN];

%one integration starting from second point
[~, output1] = ode113(@landing_dyn, [tspan_l(2) tspan_l(N)], guess(step_var+1:step_var+step_st)', options, [0;0;0;0], par);

% Check guess point
[time, output] = ode113(@CR3BP_dyn, [t1*TU/TT 24*TU/TT], state_i_rot(1:6)', options, mass_ratio); 
r_enc = zeros(length(time), 3);
for k = 1:length(time)
    % From rotating Saturn-Enceladus to IAU_Enceladus
    state_in = rot2iau_enc(time(k), output(k,:), mass_ratio);
    r_in = state_in(1:3)*DD/DU;
    v_in = state_in(4:6)*(DD/TT)/VU;
    
    % DCM Matrix: rigid rotation around X
    A_rotx = [1  0  0
              0  0  1
              0 -1  0];
    r_enc(k,:) = A_rotx*r_in;
end

Enceladus_3D(1, [0 0 0]);
hold on; grid on; grid minor
for k = 1:N
    plot3(guess((k-1)*step_var+1), guess((k-1)*step_var+2), guess((k-1)*step_var+3), '.r', 'LineWidth', 1.2);
    axis equal
end
plot3(r_enc(:,1), r_enc(:,2), r_enc(:,3), '-m', 'LineWidth', 0.5);
% plot3(output1(:,1), output1(:,2), output1(:,3), '.r', 'LineWidth', 1.2)
title('Initial Guess')
xlabel('$x$')
ylabel('$y$')

% check if initial guess satisfies constraints
for k = 1:N
    mass_check(k) = guess(step_var*(k-1)+step_st);
    u_check(k) = guess(step_var*(k-1)+(step_st+1));
end
flag_mass = find(mass_check < m_dry);
flag_umin = find(u_check < 0);
flag_umax = find(u_check > 1);


%%
% constraints for fmincon
A = [zeros(step_var*N,1); 1; -1]';
b = 0;
Aeq = [];
beq = [];

lb = -Inf*ones(1,step_var*N+2);
ub = Inf*ones(1,step_var*N+2);
for k = 1:N
    lb(step_var*(k-1)+step_st) = m_dry;
    lb(step_var*(k-1)+(step_st+1)) = 0;
    ub(step_var*(k-1)+(step_st+1)) = 1;
end
% lb(end-1) = 0;
lb(end) = 0;

%%

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp',...
    'SpecifyObjectiveGradient', false, 'MaxIter', 1000, 'MaxFunctionEvaluations', 3*1e5);

%[x_final, fval, exitflag, struct] = fmincon(@(var) land_objfun(var, state_i_rot, par, N),guess,A,b,Aeq,beq,lb,ub, ...
%    @(var) land_nonlincon(var, state_i_rot, par, N),options);
%
 %[x_final, fval, exitflag, struct] = fmincon(@(var) land_objfun(var, state_i_rot,lonlat, par, N),guess,A,b,Aeq,beq,lb,ub, ...
  %  @(var) land_nonlincon_pp(var, state_i_rot, lonlat, par, N),options);
[x_final, fval, exitflag, struct] = fmincon(@(var) land_objfunpp(var, state_i_rot,lonlat, par, N),guess,A,b,Aeq,beq,lb,ub, ...
    @(var) land_nonlincon(var, state_i_rot, par, N),options);



%% Pin-point: plot target landing site
% initial and final time
t1 = x_final(end-1);
tN = x_final(end);

% Point latitude and longitude
lat = deg2rad(lonlat(1));
lon = deg2rad(lonlat(2));
we = par(7);
Re = par(5);

% From latitudinal to cartesian
r_xz = Re*cos(lat);
Xrot = r_xz*sin(lon);
Yrot = Re*sin(lat);
Zrot = r_xz*cos(lon);
vec_rot = [Xrot; Yrot; Zrot];

% Enceladus rotation
th_e = we*(tN-t1) + th_e0;

% From rotating enceladus to IAU_Enceladus
A_rot2IAU = [cos(th_e)    0     sin(th_e)
                 0        1         0
            -sin(th_e)    0     cos(th_e)];
% desired final state
vec_pp = A_rot2IAU*vec_rot;

% final distance wrt desired landing site
rr_fin = x_final(end-12:end-10);
dist = norm(vec_pp-rr_fin);

%%
% Trajectory
Enceladus_3D(1, [0 0 0]);
hold on; grid on; grid minor
for k = 1:N
    plot3(x_final((k-1)*step_var+1), x_final((k-1)*step_var+2), x_final((k-1)*step_var+3),'.r', 'LineWidth', 2);
    plot3(r_enc(:,1), r_enc(:,2), r_enc(:,3), '-b', 'LineWidth', 2);
    %plot3(Re*exp(1i*circ),'-b')
    axis equal
end
title('Optimized Landing Trajectory')
plot3(vec_pp(1),vec_pp(2),vec_pp(3),'*')
xlabel('$x$')
ylabel('$y$')


% Control law
figure; hold on; grid on; grid minor
t_plt = linspace(x_final(end-1), x_final(end), N);
for k = 1:N
    control(k) = x_final((k-1)*step_var+step_st+1);
    plot(t_plt(k),control(k), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
end
title('Control Law')
xlabel('$t\ [hours]$')
ylabel('$u$')

% Thrust
Thrust_min = min(control)*(Tmax*FU)*1e3;     %[N]
Thrust_max = max(control)*(Tmax*FU)*1e3;     %[N]

% Final Velocity
vv_fin = x_final(end-9:end-7)*VU;
v_fin = norm(vv_fin)*1e3;                    %[m/s]
vel = zeros(N,3);
vel_norm = zeros(N,1);
for k =1:N
    vel(k,:) = x_final((k-1)*step_var+4:(k-1)*step_var+6);
    vel_norm(k) = norm(vel(k,:)); 
end
figure; hold on; grid on; grid minor
plot(t_plt, vel_norm*VU, 'ob', 'LineWidth', 1.2)
title('Velocity')
xlabel('$time\ [hours]$')
ylabel('$v\ [km/s]$')

% Mass
figure; hold on; grid on; grid minor
t_plt = linspace(x_final(end-1), x_final(end), N);
for k = 1:N
    mass(k) = x_final((k-1)*step_var+step_st);
    plot(t_plt(k),mass(k)*MM, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
end
title('Mass variation')
xlabel('$t\ [hours]$')
ylabel('$M\ [kg]$')


% Propellant Mass
m_fin = x_final(end-6)*MM;
m0_dim = m0*MM;
m_prop = m0_dim - m_fin;

% Initial thrust for first "manoeuvre"
init_thrust = 0;
for i = 1:length(control)
    if control(i) > 1e-10 
        init_thrust = control(i)*(Tmax*FU) + init_thrust;
    else
        break
    end
end


% DV: Tsiolkowsky
Is_dim = Isp*TU;
g0_dim = g0*acc;
dV = Is_dim*g0_dim*log(m0_dim/m_fin);

% Print Results
fprintf('RESULTS:\n\nMinimum thrust: %.4f [N]\nFinal velocity: %e [m/s]\nPropellant mass: %.1f [kg]\nDelta V: %f [km/s]', Thrust_min, v_fin, m_prop, dV);


%%
% validation and plot
% propagate initial orbit from t0 to t1
% propagate vector of NLP variables and check consistency

Enceladus_3D(1, [0 0 0]);
hold on; grid on; grid minor
t1 = x_final(end-1);
tN = x_final(end);
%plot3(r_enc(:,1), r_enc(:,2), r_enc(:,3), '--k', 'LineWidth', 1); % initial orbit
[time_enc, output_enc] = ode113(@CR3BP_dyn, [t1*TU/TT 24*TU/TT], state_i_rot(1:6)', options, mass_ratio); 
r_enc_plot = zeros(length(time), 3);
for k = 1:length(time_enc)
    % From rotating Saturn-Enceladus to IAU_Enceladus
    state_in = rot2iau_enc(time_enc(k), output_enc(k,:), mass_ratio);
    r_in = state_in(1:3)*DD/DU;
    v_in = state_in(4:6)*(DD/TT)/VU;
    
    % DCM Matrix: rigid rotation around X
    A_rotx = [1  0  0
              0  0  1
              0 -1  0];
    r_enc_plot(k,:) = A_rotx*r_in;
end


plot3(r_enc(:,1), -r_enc(:,3), r_enc(:,2),'LineWidth',2); % initial orbit
plot3(vec_pp(1),-vec_pp(3),vec_pp(2),'*m','LineWidth',2.5,'Markersize',8)
% [~, output_initial_orbit] = ode113(@dyn, [0 t1], [r_i;v_i], options, par);
% circ_admissible = deg2rad(180+70):pi/1e4:deg2rad(270+20);
% landing_site = Re*exp(1i*circ_admissible);
% plot(landing_site,'-g','LineWidth',10)


s0 = x_final(1:6);
view(3);

tspan = linspace(t1,tN,N);
u_plot = x_final(step_st+1);
x_plot = []; y_plot = []; z_plot = [];
for k = 1:N-1
    s0 = x_final((k-1)*step_var+1:(k-1)*step_var+step_st);
    u_plot = x_final((k-1)*step_var+step_st+1:(k-1)*step_var+step_st+4);
    [~, output] = ode113(@landing_dyn, [tspan_l(k) tspan_l(k+1)], s0, options, u_plot, par);
%     x_plot = [x_plot;output(:,1)];
%     y_plot = [y_plot;output(:,2)];
%     z_plot = [z_plot;output(:,3)];
    y_plot = [y_plot;-output(:,3)];
    z_plot = [z_plot;output(:,2)];
    x_plot = [x_plot;output(:,1)];
end
plot3(x_plot,y_plot,z_plot,'LineWidth',2,'color','r')

for k = 1:N
%     plot3(x_final((k-1)*step_var+1), x_final((k-1)*step_var+2), x_final((k-1)*step_var+3), 'om', 'LineWidth', 1.2,'MarkerSize',5);
 plot3(x_final((k-1)*step_var+1), -x_final((k-1)*step_var+3), x_final((k-1)*step_var+2), 'ok', 'LineWidth', 1.5,'MarkerSize',4,'markerfacecolor','k');
    axis equal
end
plot3(vec_pp(1),-vec_pp(3),vec_pp(2),'*m','LineWidth',2.5,'Markersize',8)

title('Fuel-optimal Landing Trajectory. $DU = 251.1\ km$')
xlabel('$x\ [DU]$')
ylabel('$y\ [DU]$')
zlabel('$z\ [DU]$')
legend('Enceladus','Initial science orbit','Target Landing Site','Landing trajectory','NLP points')

%%
% TO DO:
% Change initial guess (?)
% Insert attitude

% ASK Guadagnini:
% Downrange as path constraint
% Initial guess: Do we need to divide landing if altitude is high?


%% Try pinpoint landing as a hard constraint with different initial guess
t1 = eps;
tN = t1+0.7;
lonlat = [-70; 270]; th_e0 = deg2rad(270);
par(8) = th_e0 ;
lat = deg2rad(lonlat(1));
lon = deg2rad(lonlat(2));
% From latitudinal to cartesian
r_xz = Re*cos(lat);
Xrot = r_xz*sin(lon);
Yrot = Re*sin(lat);
Zrot = r_xz*cos(lon);
vec_rot = [Xrot; Yrot; Zrot];

% Enceladus rotation
th_e = we*(tN-t1);

% From rotating enceladus to IAU_Enceladus
A_rot2IAU = [cos(th_e)    0     sin(th_e)
                 0        1         0
            -sin(th_e)    0     cos(th_e)];
% desired final state
vec_pp = A_rot2IAU*vec_rot;
% Time grid and initial guess fill
h = (tN-t1)/(N-1);
s0 = state_i;
guess = zeros(step_var*N+2, 1);
guess(1:step_st,1) = s0;
guess(step_st+1) = 1-eps;
guess(step_st+2:step_st+4) = -s0(4:6)'/norm(s0(4:6));
guess(end-12:end-10) = vec_pp;
guess(end-5) =  1-eps;
guess(end-4:end-2) = -(s0(1:3)'-vec_pp)/norm(s0(1:3)'-vec_pp);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

tspan_l(1) = t1;
tspan_l(N) = tN;

for k = 1:N-1
    guess(k*step_var+1:k*step_var+3) = s0(1:3)'+k*h*(vec_pp-s0(1:3)');
end
guess(end-1:end) = [t1; tN];


% Check guess point
[time, output] = ode113(@CR3BP_dyn, [t1*TU/TT 24*TU/TT], state_i_rot(1:6)', options, mass_ratio); 
r_enc = zeros(length(time), 3);
for k = 1:length(time)
    % From rotating Saturn-Enceladus to IAU_Enceladus
    state_in = rot2iau_enc(time(k), output(k,:), mass_ratio);
    r_in = state_in(1:3)*DD/DU;
    v_in = state_in(4:6)*(DD/TT)/VU;
    
    % DCM Matrix: rigid rotation around X
    A_rotx = [1  0  0
              0  0  1
              0 -1  0];
    r_enc(k,:) = A_rotx*r_in;
end
figure
% Enceladus_3D(1, [0 0 0]);
hold on; grid on; grid minor
for k = 1:N
    plot3(guess((k-1)*step_var+1), guess((k-1)*step_var+2), guess((k-1)*step_var+3), '.r', 'LineWidth', 1.2);
    axis equal
end
plot3(r_enc(:,1), r_enc(:,2), r_enc(:,3), '-m', 'LineWidth', 0.5);
% plot3(output1(:,1), output1(:,2), output1(:,3), '.r', 'LineWidth', 1.2)
title('Initial Guess')
xlabel('$x$')
ylabel('$y$')

% solve
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp',...
    'SpecifyObjectiveGradient', false, 'MaxIter', 1000, 'MaxFunctionEvaluations', 3*1e5);

%[x_final, fval, exitflag, struct] = fmincon(@(var) land_objfun(var, state_i_rot, par, N),guess,A,b,Aeq,beq,lb,ub, ...
%    @(var) land_nonlincon(var, state_i_rot, par, N),options);
%
[x_final, fval, exitflag, struct] = fmincon(@(var) land_objfun(var, state_i_rot, par, N),guess,A,b,Aeq,beq,lb,ub, ...
   @(var) land_nonlincon_pp(var, state_i_rot, lonlat, par, N),options);

%% Pin-point: plot target landing site
% initial and final time
t1 = x_final(end-1);
tN = x_final(end);

% Point latitude and longitude
lat = deg2rad(lonlat(1));
lon = deg2rad(lonlat(2));
we = par(7);
Re = par(5);

% From latitudinal to cartesian
r_xz = Re*cos(lat);
Xrot = r_xz*sin(lon);
Yrot = Re*sin(lat);
Zrot = r_xz*cos(lon);
vec_rot = [Xrot; Yrot; Zrot];

% Enceladus rotation
th_e = we*(tN-t1) + th_e0;
% From rotating enceladus to IAU_Enceladus
A_rot2IAU = [cos(th_e)    0     sin(th_e)
                 0        1         0
            -sin(th_e)    0     cos(th_e)];
% desired final state
vec_pp = A_rot2IAU*vec_rot;

% final distance wrt desired landing site
rr_fin = x_final(end-12:end-10);
dist = norm(vec_pp-rr_fin);

%%
% Trajectory
Enceladus_3D(1, [0 0 0]);
hold on; grid on; grid minor
for k = 1:N
    plot3(x_final((k-1)*step_var+1), x_final((k-1)*step_var+2), x_final((k-1)*step_var+3),'.r', 'LineWidth', 2);
    plot3(r_enc(:,1), r_enc(:,2), r_enc(:,3), '-b', 'LineWidth', 2);
    %plot3(Re*exp(1i*circ),'-b')
    axis equal
end
title('Optimized Landing Trajectory')
plot3(vec_pp(1),vec_pp(2),vec_pp(3),'*')
xlabel('$x$')
ylabel('$y$')


% Control law
figure; hold on; grid on; grid minor
t_plt = linspace(x_final(end-1), x_final(end), N);
for k = 1:N
    control(k) = x_final((k-1)*step_var+step_st+1);
    plot(t_plt(k),control(k), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
end
title('Control Law')
xlabel('$t\ [hours]$')
ylabel('$u$')

% Thrust
Thrust_min = min(control)*(Tmax*FU)*1e3;     %[N]
Thrust_max = max(control)*(Tmax*FU)*1e3;     %[N]

% Final Velocity
vv_fin = x_final(end-9:end-7)*VU;
v_fin = norm(vv_fin)*1e3;                    %[m/s]
vel = zeros(N,3);
vel_norm = zeros(N,1);
for k =1:N
    vel(k,:) = x_final((k-1)*step_var+4:(k-1)*step_var+6);
    vel_norm(k) = norm(vel(k,:)); 
end
figure; hold on; grid on; grid minor
plot(t_plt, vel_norm*VU, 'ob', 'LineWidth', 1.2)
title('Velocity')
xlabel('$time\ [hours]$')
ylabel('$v\ [km/s]$')

% Mass
figure; hold on; grid on; grid minor
t_plt = linspace(x_final(end-1), x_final(end), N);
for k = 1:N
    mass(k) = x_final((k-1)*step_var+step_st);
    plot(t_plt(k),mass(k)*MM, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
end
title('Mass variation')
xlabel('$t\ [hours]$')
ylabel('$M\ [kg]$')


% Propellant Mass
m_fin = x_final(end-6)*MM;
m0_dim = m0*MM;
m_prop = m0_dim - m_fin;

% Initial thrust for first "manoeuvre"
init_thrust = 0;
for i = 1:length(control)
    if control(i) > 1e-10 
        init_thrust = control(i)*(Tmax*FU) + init_thrust;
    else
        break
    end
end


% DV: Tsiolkowsky
Is_dim = Isp*TU;
g0_dim = g0*acc;
dV = Is_dim*g0_dim*log(m0_dim/m_fin);

% Print Results
fprintf('RESULTS:\n\nMinimum thrust: %.4f [N]\nFinal velocity: %e [m/s]\nPropellant mass: %.1f [kg]\nDelta V: %f [km/s]', Thrust_min, v_fin, m_prop, dV);


%%
% validation and plot
% propagate initial orbit from t0 to t1
% propagate vector of NLP variables and check consistency

Enceladus_3D(1, [0 0 0]);
hold on; grid on; grid minor
t1 = x_final(end-1);
tN = x_final(end);
%plot3(r_enc(:,1), r_enc(:,2), r_enc(:,3), '--k', 'LineWidth', 1); % initial orbit
[time_enc, output_enc] = ode113(@CR3BP_dyn, [t1*TU/TT 24*TU/TT], state_i_rot(1:6)', options, mass_ratio); 
r_enc_plot = zeros(length(time), 3);
for k = 1:length(time_enc)
    % From rotating Saturn-Enceladus to IAU_Enceladus
    state_in = rot2iau_enc(time_enc(k), output_enc(k,:), mass_ratio);
    r_in = state_in(1:3)*DD/DU;
    v_in = state_in(4:6)*(DD/TT)/VU;
    
    % DCM Matrix: rigid rotation around X
    A_rotx = [1  0  0
              0  0  1
              0 -1  0];
    r_enc_plot(k,:) = A_rotx*r_in;
end


plot3(r_enc(:,1), -r_enc(:,3), r_enc(:,2),'LineWidth',2); % initial orbit
plot3(vec_pp(1),-vec_pp(3),vec_pp(2),'*m','LineWidth',2.5,'Markersize',8)
% [~, output_initial_orbit] = ode113(@dyn, [0 t1], [r_i;v_i], options, par);
% circ_admissible = deg2rad(180+70):pi/1e4:deg2rad(270+20);
% landing_site = Re*exp(1i*circ_admissible);
% plot(landing_site,'-g','LineWidth',10)


s0 = x_final(1:6);
view(3);

tspan = linspace(t1,tN,N);
u_plot = x_final(step_st+1);
x_plot = []; y_plot = []; z_plot = [];
for k = 1:N-1
    s0 = x_final((k-1)*step_var+1:(k-1)*step_var+step_st);
    u_plot = x_final((k-1)*step_var+step_st+1:(k-1)*step_var+step_st+4);
    [~, output] = ode113(@landing_dyn, [tspan_l(k) tspan_l(k+1)], s0, options, u_plot, par);
%     x_plot = [x_plot;output(:,1)];
%     y_plot = [y_plot;output(:,2)];
%     z_plot = [z_plot;output(:,3)];
    y_plot = [y_plot;-output(:,3)];
    z_plot = [z_plot;output(:,2)];
    x_plot = [x_plot;output(:,1)];
end
plot3(x_plot,y_plot,z_plot,'LineWidth',2,'color','r')

for k = 1:N
%     plot3(x_final((k-1)*step_var+1), x_final((k-1)*step_var+2), x_final((k-1)*step_var+3), 'om', 'LineWidth', 1.2,'MarkerSize',5);
 plot3(x_final((k-1)*step_var+1), -x_final((k-1)*step_var+3), x_final((k-1)*step_var+2), 'ok', 'LineWidth', 1.5,'MarkerSize',4,'markerfacecolor','k');
    axis equal
end
plot3(vec_pp(1),-vec_pp(3),vec_pp(2),'*m','LineWidth',2.5,'Markersize',8)

title('Fuel-optimal Landing Trajectory. $DU = 251.1\ km$')
xlabel('$x\ [DU]$')
ylabel('$y\ [DU]$')
zlabel('$z\ [DU]$')
legend('Enceladus','Initial science orbit','Target Landing Site','Trajectory','NLP points')
 