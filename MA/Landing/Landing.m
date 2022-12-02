%% Landing problem: Collocation method
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(0,'defaultAxesFontSize', 16)
%

%%
clearvars; close all; 
clc 

% step1: vehicle par (TO CHANGE!!!)
Tmax = 4*125e-3;                %[kN] Maximum Thrust  !!!!! 4* %%%%%%%%%%%%%%%%
Isp = 228;                      %[s] Specific Impulse
g0 = 9.81*1e-3;                 %[km/s^2] acceleration constant
% m0 = 95+75.4;                   %[kg] initial mass of lander (both Non Sampling Orbiter- Sampling Lander and S-S)
% m_dry = 75.4;

m0 = 95+746;                   % [kg] initial mass of lander ( Non Sampling Orbiter- Sampling Lander)
m_dry = 746;


% Enceladus par
Re = 251.1;                                       %[km] mean radius of Enceladus
mu = (6.67430e-11 *8.6e19 )*10^(-9) ;             % [km^3/s^2] Enceladus gravitational constant
mass_ratio = 1.90095713928102*1e-7;               % Saturn-Enceladus mass ratio

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

% organize par 
par(1) = Tmax;
par(2) = Isp;
par(3) = g0;
par(4) = mu;
par(5) = Re;
par(6) = mass_ratio;

N = 50;

% INITIAL ORBIT : circular, polar (TO CHANGE!!)
initial = 4;
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
        h = 100/DU;
        r_mod = Re+h;
        xi = -r_mod;
        yi = 0;
        %zi = ;
        v_mod = sqrt(mu/norm([xi; yi]));
        vxi = 0;
        vyi = -v_mod;
        %vzi = ;
        %state_i = [xi yi zi vxi vyi vzi m0];
        % Initial time guess
        t1 = eps; 
        % Final time guess if h = 100 km
        tN = 0.9*3600/TU;             
    case 3
        DD=238411468.296/1000; %km
        TT=118760.57/(2*pi);
        xi = 1.000062853735440;
        yi = 0;
        zi = -0.00117884381145460;
        vxi = 0;
        vyi = 0.0168877463349484;
        vzi = 0;
        state_i_rot = [xi yi zi vxi vyi vzi m0];
        % Initial time guess
        t1 = 7*3600/TU; 
        % Final time guess 
        tN = 0.7*3600/TU + t1;   

        % From rotating Saturn-Enceladus to IAU_Enceladus
        state_in = rot2iau_enc(t1, state_i_rot(1:end-1), mass_ratio);
        r_in = state_in(1:3)*DD/DU;
        v_in = state_in(4:6)*(DD/TT)/VU;
        
        % DCM Matrix: rigid rotation around X
        A_rotx = [1  0  0
                  0  0  1
                  0 -1  0];
        r1 = A_rotx*r_in;
        v1 = A_rotx*v_in;
        state_i = [r1' v1' m0];
    case 4
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
        state_in = rot2iau_enc(t1, state_i_rot(1:end-1), mass_ratio);
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
lb(end-1) = 0;
lb(end) = 0;

%%

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp',...
    'SpecifyObjectiveGradient', false, 'MaxIter', 1000, 'MaxFunctionEvaluations', 3*1e5);

[x_final, fval, exitflag, struct] = fmincon(@(var) land_objfun(var, state_i_rot, par, N),guess,A,b,Aeq,beq,lb,ub, ...
    @(var) land_nonlincon(var, state_i_rot, par, N),options);

%%
% Trajectory
Enceladus_3D(1, [0 0 0]);
hold on; grid on; grid minor
for k = 1:N
    plot3(x_final((k-1)*step_var+1), x_final((k-1)*step_var+2), x_final((k-1)*step_var+3),'.r', 'LineWidth', 1.2);
    plot3(r_enc(:,1), r_enc(:,2), r_enc(:,3), '-m', 'LineWidth', 0.5);
    %plot3(Re*exp(1i*circ),'-b')
    axis equal
end
title('Optimized Landing Trajectory')
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
plot3(r_enc(:,1), r_enc(:,2), r_enc(:,3), '--k', 'LineWidth', 1); % initial orbit
[~, output_initial_orbit] = ode113(@dyn, [0 t1], [r_i;v_i], options, par);
circ_admissible = deg2rad(180+70):pi/1e4:deg2rad(270+20);
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
    x_plot = [x_plot;output(:,1)];
    y_plot = [y_plot;output(:,2)];
    z_plot = [z_plot;output(:,3)];
end
plot3(x_plot,y_plot,z_plot,'b','LineWidth',1.5)

for k = 1:N
    plot3(x_final((k-1)*step_var+1), x_final((k-1)*step_var+2), x_final((k-1)*step_var+3), 'om', 'LineWidth', 1.2,'MarkerSize',5);
    axis equal
end


title('Fuel-optimal Landing Trajectory. $DU = 251.1\ km$')
xlabel('$x\ [DU]$')
ylabel('$y\ [DU]$')
legend('Enceladus','Initial orbit $h = 60\ km$','Admissible landing','Landing trajectory','NLP points')

%%
% TO DO:
% Change initial guess (?)
% Insert attitude

% ASK Guadagnini:
% Downrange as path constraint
% Initial guess: Do we need to divide landing if altitude is high?



