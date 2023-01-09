%% Uncertainties propagation
clearvars; close all; clc

% 1- 50 points
% 2 - 100 points
points = 1; 
switch points
    case 1
        load 'optimal50.mat'
        N = 50;
    case 2
        load 'optimal100.mat'
        N = 100;
end

n_T = 4; % number of thrusters.
Tmax = n_T*220e-3; %[kN] Maximum Thrust 
Isp = 229; %[s] Specific Impulse
g0 = 9.81*1e-3; %[km/s^2] acceleration constant

m_dry_lander =508; % [kg]
m_prop_HA = 80; % [kg] mass of propellant for hazard avoidance. NOTE: HYP
m_dry = m_dry_lander + m_prop_HA; % NOTE: improper use of name "dry mass": 
m_prop = 102.2; % [kg] NOTE: estimated mass of propellant from PS, can be iterated
m0 = m_dry+ m_prop; % [kg] WET MASS: initial mass of landing trajectory

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

step_var = 11;
step_st = 7;

% Uncertainties: spring + sensors
x_hat0 = x_final(1:6);
% Covariance: adimensionalize with [km^2, km^2/s, km^2/s^2]
% For now: 10m and 1cm/s
P0 = [1e-4/DU^2   0            0           0               0               0 
      0           1e-4/DU^2    0           0               0               0 
      0           0            1e-4/DU^2   0               0               0 
      0           0            0           1e-10/(DU/TU)^2 0               0 
      0           0            0           0               1e-10/(DU/TU)^2 0
      0           0            0           0               0               1e-10/(DU/TU)^2];
% P0 = zeros(6,6);
t1 = x_final(end-1);
tN = x_final(end);
tspan_new = linspace(t1, tN, N);
options_ode = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

% Find 2km altitude position
for k = 1:N-1
    x_des = x_final((k-1)*step_var+1:(k-1)*step_var+3);
    alt = norm(x_des - Re.*x_des/norm(x_des));
    if alt <= 2/DU
        index_alt = k;
        break
    end
end

%% Unscented Transform
% alpha = 1e-3;
% beta = 2;
% n = 6;
% lambda = alpha^2*(n)-n;
% M = (n+lambda)*P0;
% R = sqrtm(M);
% 
% chi = zeros(6,2*n+1);  %each column is a sigma point
% Wm = zeros(1,2*n+1);
% Wc = zeros(1,2*n+1);
% for i = 1:(2*n+1)
%     if i <= n
%         chi(:,i) = x_hat0 + R(:,i);
%         Wm(i) = 1/(2*(n+lambda));
%         Wc(i) = 1/(2*(n+lambda));
%     elseif i == n+1
%         chi(:,i) = x_hat0;
%         Wm(i) = lambda/(n+lambda);
%         Wc(i) = lambda/(n+lambda) + (1-alpha^2 + beta);
%     else
%         chi(:,i) = x_hat0 - R(:,i-(n+1));
%         Wm(i) = 1/(2*(n+lambda));
%         Wc(i) = 1/(2*(n+lambda));
%     end
% end
% chi(end+1,:) = m0;
% 
% mean_ut = x_hat0;
% P_ut = P0;
% y = zeros(6,2*n+1);
% for j = 1:(2*n+1)
%     state0 = chi(:,j);
%     for k = 1:index_alt-1
%         u_plot = x_final((k-1)*step_var+step_st+1:(k-1)*step_var+step_st+4);
%         [~, output] = ode113(@landing_dyn, [tspan_new(k) tspan_new(k+1)], state0, options_ode, u_plot, par);
%         state0 = output(end,1:7)';
%     end
%     y(:,j) = state0(1:6,1);
%     mean_ut = Wm(j).*y(:,j) + mean_ut;
% end  
% for j = 1:(2*n+1)
%     col = y(:,j) - mean_ut;
%     row = col';
%     P_ut = (Wc(j) .* (col*row)) + P_ut;
% end

%% Montecarlo 
N_mc = 500;
P_mc = P0;
for k = 1:index_alt-1 
    x_hat = x_final((k-1)*step_var+1:(k-1)*step_var+6);
    samples = mvnrnd(x_hat, P_mc, N_mc);
    state_mc = zeros(N_mc,6);
    state0 = zeros(1,7);
    for i=1:N_mc
        state0(1:6) = samples(i,:);
        state0(end) = x_final((k-1)*step_var+7);
        state0 = state0';
        u_plot = x_final((k-1)*step_var+step_st+1:(k-1)*step_var+step_st+4);
        [~, output] = ode113(@landing_dyn, [tspan_new(k) tspan_new(k+1)], state0, options_ode, u_plot, par);
        state0 = output(end,1:7)';
        state_mc(i,:) = state0(1:6);
    end
    mean_mc = mean(state_mc, 1);
    P_mc = cov(state_mc);

end

% Plot ellipse 
p = 0.997;
Enceladus_3D(1, [0 0 0])
% figure;
% e1 = plot_gaussian_ellipsoid(mean_ut(1:3), P_ut(1:3,1:3), 3);
hold on
e2 = plot_gaussian_ellipsoid(mean_mc(1:3), P_mc(1:3,1:3), 3);
% m1 = plot3(mean_ut(1), mean_ut(2), mean_ut(3), 's', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120', 'MarkerSize', 10);
m2 = plot3(mean_mc(1), mean_mc(2), mean_mc(3), 's', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 10);
p1 = plot3(x_final((index_alt-1)*step_var+1), x_final((index_alt-1)*step_var+2), x_final((index_alt-1)*step_var+3), 'o', 'MarkerSize', 10);
s_mc = plot3(state_mc(:,1), state_mc(:,2), state_mc(:,3), '.', 'color', '#D95319');
uistack(e2, 'bottom');
xlabel('X [-]')
ylabel('Y [-]')
axis equal
s = legend([p1 e2 m2 s_mc], 'Landing Site', '$3\sigma - MC$', '$MC - \hat{r}_{xy}$', '$MC samples$');
s.FontSize = 12;
grid on; grid minor;