%% EX 1
% 3D Saturn-Enceladus CRTBP with mu = 1.513e-7;

clearvars; clc

% initial conditions
x0 = 0.997637174347;
y0 = 0;
z0 = 0.002523536751;
vx0 = 0;
vy0 = 0.002967906138;
vz0 = 0;
state0 = [x0; y0; z0; vx0; vy0; vz0];
state0_c = state0;

mu = 0.0000001480;

syms x y z vx vy vz
dxdt = [vx
        vy
        vz
        2*vy + x - (1-mu)*(x+mu)/(norm([x+mu; y; z]))^3 - mu*(x+mu-1)/(norm([x+mu-1; y; z]))^3
        -2*vx + y - (1-mu)*y/(norm([x+mu; y; z]))^3 - mu*y/(norm([x+mu-1; y; z]))^3
        -(1-mu)*z/(norm([x+mu; y; z]))^3 - mu*z/(norm([x+mu-1; y; z]))^3];

% matrix A = dfdt
v = [x; y; z; vx; vy; vz];
dfdt_sym = jacobian(dxdt, v);
dfdt = matlabFunction(dfdt_sym);

%initial STM
psi0 = eye(6,6);
psi0_shaped = reshape(psi0, 1, length(psi0(1,:))*length(psi0(:,1))); 

%integration
N = 1e4;
t = linspace(0, pi, N);
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @ObjEvent);

% initial orbit
tspan = linspace(0, pi, 100);
[time, output1] = ode113( @(tspan,s) CRTBP_integrator(tspan, s, mu, dfdt), tspan, [state0; psi0_shaped'], options);
state1 = output1(:,1:6);

% family of orbits
zf = z0;
n = 1;
zv = linspace(z0, zf, n);
err = 1;
% output_e = zeros(n,42);
% output_e(1,1:6) = state0;
tic
count = 0;
tevent = zeros(1,length(zv));
for i = 1:length(zv)
    vx_e = 1;
    vz_e = 1;
    delta_vy = 0;
    delta_x = 0;
    figure;
    while abs(vx_e) > 1e-6 || abs(vz_e) > 1e-6
        state0(1) = state0(1) + delta_x;
        state0(5) = state0(5) + delta_vy;
        [time, output, te, output_e, ie] = ode113( @(t,s) CRTBP_integrator(t, s, mu, dfdt), t, [state0; psi0_shaped'], options);
        
        vx_e = output_e(4);
        vz_e = output_e(6);
        psi_e = output_e(7:end);
        psi_e = reshape(psi_e, length(psi0(1,:)), length(psi0(:,1)));
        ratio = psi_e(6,5)/psi_e(4,5);
        M = [psi_e(4,1) psi_e(4,5)
             psi_e(6,1) psi_e(6,5)];
        v = [-vx_e; -vz_e];
        u = M\v;
        delta_x = u(1);
        delta_vy = u(2);
        count = count + 1;
        
        %plot temp
        state_plt = output(:,1:6);
        plot3(state_plt(:,1), state_plt(:,2), state_plt(:,3));
        hold on
        grid on
    end
    tevent(i) = te;
    s(i,:) = state0;
    if i ~= length(zv)
        state0(3) = zv(i+1);
    end
end
fprintf('Computational Time: %.2f', toc);
fprintf('\nNumber of iterations: %d', count);

%% plot orbit
close all

s = state0;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
for i=1:10:length(s(:,1))
%     tspan = linspace(0, 2*tevent(i), 1000);
    tspan = linspace(0, 2*pi, 1000);
    s1 = s(i,:)';
    s1 = s;
    [time, output_long] = ode113( @(t,s) CRTBP_integrator(tspan, s, mu, dfdt), tspan, [s1; psi0_shaped'], options);
    state_long = output_long(:,1:6);
    plot3(state_long(:,1), state_long(:,2), state_long(:,3), 'b-');
    hold on
end
% plot3(state1(:,1), state1(:,2), state1(:,3), 'r-');
title('Halo Orbits');
xlabel('x');
ylabel('y');
zlabel('z');
grid on
