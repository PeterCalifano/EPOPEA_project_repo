clear;
close all;
clc;

%% INERTIA MATRIX: Stored
% Inertia matrix of S/C with deployed solar panels
mass = 150;

%solar panels
lx_sp = 1.5;
ly_sp = 0.02;
lz_sp = 1;
dens_sp = 3.8; %[kg/m^2]
m_sp = dens_sp*ly_sp*lz_sp;

Jy_sp = (1/12)*m_sp*(lz_sp^2 + lx_sp^2);
Jz_sp = (1/12)*m_sp*(lx_sp^2 + ly_sp^2);
Jx_sp = (1/12)*m_sp*(lz_sp^2 + ly_sp^2);

%Main body
lx_mb = 1; %[m]
ly_mb = 1;
lz_mb = 1;
m_mb = mass - m_sp;

Jx_mb = (1/12)*m_mb*(ly_mb^2 + lz_mb^2);
Jy_mb = (1/12)*m_mb*(lx_mb^2 + lz_mb^2);
Jz_mb = (1/12)*m_mb*(ly_mb^2 + lx_mb^2);

%Huygens-Steiner
d = ly_mb/2 + ly_sp/2;
Jx_sp = Jx_sp + m_sp*d^2;
Jz_sp = Jz_sp + m_sp*d^2;

%INERTIA MATRIX: STORED CONF
Jx = 2*Jx_sp + Jx_mb;
Jy = 2*Jy_sp + Jy_mb;
Jz = 2*Jz_sp + Jz_mb;

J = [Jx,0,0;
     0,Jy,0;
     0,0,Jz];
 
Jinv = inv(J);


%% Main body normals
n1 = [1;0;0]; n4 = -n1; % x direction
n2 = [0;1;0]; n5 = -n2; % y direction
n3 = [0;0;1]; n6 = -n3; % z direction

% Main Body properties
A_cube = 1;
rho_s = 0.5;
rho_d = 0.1;
r_cube = 0.5; % only valid for a 1x1x1 cube
r_i = r_cube * [n1,n2,n3,n4,n5,n6];

% Solar panels normals
n7 = [0;0;0]; n8 = -n7;
n9 = [0;0;0]; n10 = -n9;

% Solar panels CM direction
r7 = [0;0;0]; r8 = -r7;
r9 = [0;0;0]; r10 = r9;

% Solar Panel properties
A_sp = 1.5;
rho_s_sp = 0.1;
rho_d_sp = 0.1;
r_sp = ly_mb/2 + ly_sp/2;
r_i_sp = r_sp * [n7,n8,n9,n10];

%% SRP Torque
n_sun = 2*pi/(365.25*24*3600);

eps = deg2rad(23.45);
ceps = cos(eps);
seps = sin(eps);

F_e = astroConstants(31); %W/m^2
c = 299792458; %m/s
P = F_e/c; %Pa


%% Initial Attitude Data
mu = astroConstants(13);
R_earth = astroConstants(23);
a = R_earth+1200; 
e = 0; 
i = 86.4;

% Initial orbit conditions
[r0,v0] = kep2car([a,e,i,0,0,0],mu);
T = 2*pi*sqrt( a^3/mu ); % Orbital period [s]
n = 2*pi/T;

w0 = [0.3;0.2;0.5];

%% DCM Conditions
A = eye(3);
q1 = 0.5*(1+A(1,1)-A(2,2)-A(3,3))^0.5;
q2 = 0.5*(1-A(1,1)+A(2,2)-A(3,3))^0.5;
q3 = 0.5*(1-A(1,1)-A(2,2)+A(3,3))^0.5;
q4 = 0.5*(1+A(1,1)+A(2,2)+A(3,3))^0.5;
qmax = max([q1,q2,q3,q4]);

if qmax == q1
    q = [A(1,2)+A(2,1);
         A(1,3)+A(3,1);
         A(2,3)-A(3,2)]/(4*q1);
    q0 = [q1; q];
end

if qmax == q2 
   q0 = [A(1,2)+A(2,1);
         (q2^2)*4;
          A(2,3)+A(3,2);
          A(3,1)-A(1,3)]/(4*q2);
end

if qmax == q3
    q0 = [A(1,3)+A(3,1);
          A(2,3)+A(3,2);
          (q3^2)*4;
          A(1,2)-A(2,1)]/(4*q3);
end

if qmax == q4
    q = [A(2,3)-A(3,2);
         A(3,1)-A(1,3);
         A(1,2)-A(2,1)]/(4*q4);
    q0 = [q; q4];
end   


%% Magnetic Torque
% From data
m = [0.1; 0.1; 0.1]; %Am^2

% to define m_hat matrix
m_i = deg2rad(11.5); % magnetic poles inclination
w_earth = 2*pi/(24*3600);
s_m_i = sin(m_i);
c_m_i = cos(m_i);

% from IGRF 2000 (given as X_nm)
g_10 = -29615;
g_11 = -1728;
h_11 = 5186;
H0 = sqrt(g_10^2+g_11^2+h_11^2);


%% Star sensor
% Star Tracker Extended NST
% attitude error: 6 arcsec (40 arcsec about boresight)
% field of view: 10° x 12°
% mass: 1.3 kg

% Random star directions:
sy = cos(10*pi/180);
sx = sin(10*pi/180);
sx2= sin(3*pi/180);
sz = sin(5*pi/180);
sz2= sin(8*pi/180);

S_star1 = [sx;sy;0];
S_star1 = S_star1/norm(S_star1);
S_star2 = [sx2;sy;0];
S_star2 = S_star2/norm(S_star2);
S_star3 = [0;sy;sz];
S_star3 = S_star3/norm(S_star3);
S_star4 = [0;sy;-sz2];
S_star4 = S_star4/norm(S_star4);

% B* matrix
B = [S_star1,S_star2,S_star3,S_star4];
B_ast = pinv(B);

% Error matrix
alpha_x = deg2rad(6/3600); 
alpha_y = deg2rad(6/3600);
alpha_z = deg2rad(40/3600);
   
% Sensor orientation in body frame  
A_S2B = [1,0,0;0,-1,0;0,0,-1];
 

%% Actuators
% 8 actuators cofiguration, on the upper face of the cube (positive z direction) attached in pairs to th
% edges (at a distance 2d from each other, symmetrically with respect to
% the center line). Looking from the upper face of the cube (z pointing
% outwards) the pairs are numbered 1-2, 3-4, 5-6, 7-8 starting from the one pointed in 
% the positive x direction & in the positive xy quadrant and going clockwise 

% 1N Arianespace Hydrazine thruster, MIB=0.01 to 0.043 Ns
F_min = 0.160;  % N
F_max = 0.650; % N
MIB_real=F_min*0.1;

MIB_max=0.043;
Dt_min=MIB_max/F_min;

R = diag([lx_mb/2 ly_mb/2 lz_mb/2]) * [0  0  1  1  0  0 -1 -1;
                                       1  1  0  0 -1 -1  0  0;
                                      -1  1 -1  1 -1  1 -1  1];

% alpha = deg2rad(30);
% ca = cos(alpha);
% sa = sin(alpha);
% R = lx_mb/2*diag([ca-sa sa ca]) * [1 -1 1 -1;
%                                    -1 1 1 -1;
%                                    -1 -1 1 1];

w = ones(length(R(1,:)),1);

R_ast = pinv(R);


%% Control - detumbling
k_p = 0;
k_d = 100;

w_c = [0 0 0]';


%% Post Processing
% ErrP = out.ErrorAngle.Data;
% omega = out.w_B;
% acc = out.w_B_dot.Data;
% tout = out.tout;
% 
% ErrP_m = mean(ErrP)*ones(length(tout),1);
% ErrP_v = var(ErrP)*ones(length(tout),1);
% 
% figure(1);
% plot(tout,ErrP);
% hold on;
% plot(tout,ErrP_m,'k--');
% xlabel('Time [s]');
% ylabel('Error Angle [deg]');
% title('Pointing Error');
% grid on;
% 
% figure(2);
% plot(tout,acc(:,1),tout,acc(:,2),tout,acc(:,3), 'linewidth', 1);
% xlabel('Time [s]');
% ylabel('Angular acceleration [rad/s^2]');
% title('Angular Acceleration in Body Frame');
% legend('a_x','a_y','a_z');
% grid on;
% axis([0 60 -0.03 0.005]);
% 
% figure(3);
% plot(tout,omega(:,1),tout,omega(:,2),tout,omega(:,3),'linewidth', 1);
% xlabel('Time [s]');
% ylabel('Angular velocity [rad/s]');
% title('Angular Velocity in Body Frame');
% legend('\omega_x','\omega_y','\omega_z');
% grid on;
% axis([0 60 -0.1 0.6]);
