%% Relative Position
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(0,'defaultAxesFontSize', 16)

clearvars; clc; close all; cspice_kclear();

% kernelpool = fullfile('EPOPEA_metakernel.tm'); % for Windows
% cspice_furnsh(kernelpool);

% Define Kernel
cd('..')
cspice_furnsh('spice_kernels/naif0012.tls')
cspice_furnsh('spice_kernels/de440s.bsp')
cspice_furnsh('spice_kernels/sat441.bsp')
cspice_furnsh('spice_kernels/pck00010.tpc')
cspice_furnsh('spice_kernels/gm_de431.tpc')
cspice_furnsh('spice_kernels/dss.bsp')
cspice_furnsh('spice_kernels/dss.tf')
cspice_furnsh('spice_kernels/plu058.bsp')
cspice_furnsh ('spice_kernels\LATs.bsp');
cspice_furnsh ('spice_kernels\LATs.tf');
cd('MA')

%% DATA
mu=1.90095713928102*1e-7; % Saturn-Enceladus 3BP constant
DU=238411468.296/1000; %km
TU=118760.57/(2*pi); 

% Define useful constants
R_Enc = mean(cspice_bodvrd('602','RADII',3));
R_Sat = mean(cspice_bodvrd('699','RADII',3));
mu_Sun=cspice_bodvrd('SUN','GM',1);

%% FIND INITIAL TIME OF PROPAGATION

% Select the initial time to look for the instant of intersection
t0 = cspice_str2et('2051 DEC 25 12:00:00.00')/TU;
fprintf('\nInitial propagation date:\n%s', cspice_et2utc(t0*TU, 'C', 0));

% Equatorial rotation
i_ax_Sat=(-26.73)*pi/180; 

% Inclination Enceladus
i_EncSat=deg2rad(-0.009);
R_i_EncSat=[1 0 0;
            0 cos(i_EncSat) sin(i_EncSat);
            0 -sin(i_EncSat) cos(i_EncSat)];

% DEFINE the n of days to propagate
n_days = 1;

% DEFINE the number of points
n_points = 2000;

% Time Grid
tf = t0 + n_days*24*3600/TU;
tt = linspace(t0,tf,n_points);

% Initial state for the Halo orbit - Pericenter
x0_Halo=1.000062853735440;
y0_Halo=0;
z0_Halo=-0.00117884381145460;
vx0_Halo=0;
vy0_Halo=0.0168877463349484;
vz0_Halo=0;
state0_Halo=[x0_Halo,y0_Halo,z0_Halo,vx0_Halo,vy0_Halo,vz0_Halo]';

% Propagate Halo
options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);
[t_vec_Halo,state_vec_Halo]=ode113(@(t,x) CR3BP_dyn(t,x,mu),tt,state0_Halo,options_ode);
state_vec_Halo=state_vec_Halo';
x_CR3BP = state_vec_Halo(1:3,:);

% Initialize storage variables
x_sc_IS = zeros(3,length(tt));
Sat2Enc_IS = zeros(3,length(tt));
Sat2Enc_real_EC = zeros(3,length(tt));
Sat2Enc_real_IS = zeros(3,length(tt));
Sun_dir_EC = zeros(3,length(tt));
Sun_dir_IS = zeros(3,length(tt));
sc2enc = zeros(3,length(tt));
sc2enc_dir = zeros(3,length(tt));
sc_tang_dir = zeros(3,length(tt));
i_nad = zeros(1,length(tt));
i_tan = zeros(1,length(tt));
for j = 1:length(tt)

    % Real enceladus ephemerides - ECLIPTIC
    Sat2Enc_real_EC(:,j) = cspice_spkpos('602', tt(j)*TU, 'ECLIPJ2000', 'NONE', '699');

    % Saturn orbital parameters recovery
    Sun2Sat = cspice_spkezr('SATURN',  tt(j)*TU, 'ECLIPJ2000', 'NONE', 'SUN');
    r_Sat = Sun2Sat(1:3);
    v_Sat = Sun2Sat(4:6);
    [a_Sat,e_Sat,i_Sat,OM_Sat,om_Sat,theta_Sat] = car2kep(r_Sat,v_Sat,mu_Sun);
    
    % Sun wrt Saturn - ECLIPTIC
    Sat2Sun = cspice_spkpos('SUN',  tt(j)*TU, 'ECLIPJ2000', 'NONE', '699');
    Sun_dir_EC(:,j) = Sat2Sun/norm(Sat2Sun);

    % Rotation from ECLIPTIC to SATURN EQUATORIAL   
    %z-Rotation, om+theta
    R_an_Sat=[cos(om_Sat+theta_Sat) sin(om_Sat+theta_Sat) 0
             -sin(om_Sat+theta_Sat) cos(om_Sat+theta_Sat) 0
              0                     0                     1];
    %x-Rotation, i
    R_i_Sat=[1 0 0
             0 cos(i_Sat) sin(i_Sat)
             0 -sin(i_Sat) cos(i_Sat)];
    %z-Rotation OM
    R_OM_Sat=[cos(OM_Sat) sin(OM_Sat) 0
             -sin(OM_Sat) cos(OM_Sat) 0
              0 0 1];

    % x equator - Rotation 
    R_eq_Sat=[1 0 0
              0 cos(i_ax_Sat) sin(i_ax_Sat)       
              0 -sin(i_ax_Sat) cos(i_ax_Sat)];
    
    % Enceladus position in SATURN EQUATORIAL INERTIAL
    x_sc_IS(1,j) = ((x_CR3BP(1,j)+mu)*cos(tt(j)-t0) - x_CR3BP(2,j)*sin(tt(j)-t0))*DU;
    x_sc_IS(2,j) = ((x_CR3BP(1,j)+mu)*sin(tt(j)-t0) + x_CR3BP(2,j)*cos(tt(j)-t0))*DU;
    x_sc_IS(3,j) = x_CR3BP(3,j)*DU;

    % Do the same with Enceladus orbit to have a coherent reference 
    Sat2Enc_IS(1,j) = 1*cos(tt(j)-t0)*DU;
    Sat2Enc_IS(2,j) = 1*sin(tt(j)-t0)*DU;
    Sat2Enc_IS(3,j) = 0*DU;
    
    % Rotate from ECLIPTIC to INERTIAL SATURN
    Sat2Enc_real_IS(:,j) = R_i_EncSat*R_an_Sat*R_i_Sat*R_OM_Sat*R_eq_Sat*Sat2Enc_real_EC(:,j);
    Sun_dir_IS(:,j) = R_i_EncSat*R_an_Sat*R_i_Sat*R_OM_Sat*R_eq_Sat*Sun_dir_EC(:,j);

    % Enceladus wrt S/C (Nadir) - Inertial Saturn
    sc2enc(:,j) = Sat2Enc_IS(:,j) - x_sc_IS(:,j);

    % S/C (roll) - Inertial Saturn
    if sc2enc(3,j) < 0 
        sc_tang = cross(x_sc_IS(:,j),sc2enc(:,j));
    else
        sc_tang = -cross(x_sc_IS(:,j),sc2enc(:,j));
    end
    sc_tang_dir(:,j) = sc_tang/norm(sc_tang);
    sc2enc_dir(:,j) = sc2enc(:,j)/norm(sc2enc(:,j));

    %%%% TO DO %%%%
    % S/C (out-of-plane) - Inertial Saturn
    %%%% TO DO %%%%

    % RELATIVE POSITIONS OF SUN wrt NADIR, TANGENTIAL, OUT-OF-PLANE
    i_nad(j) = rad2deg(acos(dot(Sun_dir_IS(:,j), sc2enc_dir(:,j))));
    i_tan(j) = rad2deg(acos(dot(Sun_dir_IS(:,j), sc_tang_dir(:,j))));

end

% Out of plane direction
fourth = floor(length(Sat2Enc_real_IS)/4);
hh_real = cross(Sat2Enc_real_IS(:,1),Sat2Enc_real_IS(:,fourth));
h_real = hh_real/norm(hh_real);
hh_fake = cross(Sat2Enc_IS(:,1),Sat2Enc_IS(:,fourth));
h_fake = hh_fake/norm(hh_fake);

% Check angle misalignment
i_err = acos(dot(h_fake,h_real));
fprintf('\n\nOut-of-plane misalignment:\n %.5f', i_err)

% Plot without axis equal just for visualization
figure
scatter3(0,0,0,50);
hold on
plot3(Sat2Enc_IS(1,:),Sat2Enc_IS(2,:),Sat2Enc_IS(3,:),'r', 'LineWidth', 1.3)   %SPICE
plot3(x_sc_IS(1,:),x_sc_IS(2,:),x_sc_IS(3,:),'-k')                    %S/C
quiver3(zeros(1,length(tt)),zeros(1,length(tt)),zeros(1,length(tt)), ...
    1e4*Sun_dir_IS(1,:), 1e4*Sun_dir_IS(2,:), 1e4*Sun_dir_IS(3,:))
for i = 1:30:length(tt)
    quiver3(0,0,0,x_sc_IS(1,i),x_sc_IS(2,i),x_sc_IS(3,i), 'm')
    quiver3(x_sc_IS(1,i),x_sc_IS(2,i),x_sc_IS(3,i),sc2enc(1,i), sc2enc(2,i), sc2enc(3,i), 'k')
    quiver3(x_sc_IS(1,i),x_sc_IS(2,i),x_sc_IS(3,i),5e3*sc_tang_dir(1,i), 5e3*sc_tang_dir(2,i), 5e3*sc_tang_dir(3,i), 'b')
end
title(sprintf(cspice_et2utc(t0*TU, 'C', 0)))
axis equal
subtitle('Merry XMas Enceladians!')
legend('Saturn','Enceladus','S/C', 'Nadir')
xlabel('X')
ylabel('Y')
zlabel('Z')

% Angles wrt Sun
time = (tt.*TU)/86400;
shift = date2jd([2000,01,02,11,59,00])-date2jd([0000,01,01,00,00,00]);
time = time + shift;
figure;
plot(time, i_nad, 'b', 'Linewidth', 1.2)
ylabel('$\theta_{Sun/Nadir}\ [deg]$')
title('Angle between Sun and Nadir')
grid on; grid minor
ax = gca;
ax.XTick = time(1):1:time(end);
datetick('x', 'yy mmm dd - HH:MM', 'keeplimits', 'keepticks');
xtickangle(45);

figure;
plot(time, i_tan, 'r', 'Linewidth', 1.2)
ylabel('$\theta_{Sun/Tan}\ [deg]$')
title('Angle between Sun and Tangential direction')
grid on; grid minor
ax = gca;
ax.XTick = time(1):1:time(end);
datetick('x', 'yy mmm dd - HH:MM', 'keeplimits', 'keepticks');
xtickangle(45);

% % Plot for comparison
% figure
% scatter3(0,0,0,50);
% hold on
% plot3(Sat2Enc_IS(1,:),Sat2Enc_IS(2,:),Sat2Enc_IS(3,:),'r', 'LineWidth', 1.3)   %SPICE
% plot3(Sat2Enc_real_IS(1,:),Sat2Enc_real_IS(2,:),Sat2Enc_real_IS(3,:),'g--', 'LineWidth', 1.3) %REAL
% plot3(x_sc_IS(1,:),x_sc_IS(2,:),x_sc_IS(3,:),'-k')                    %S/C
% % quiver3(0,0,0,1e5*h_fake(1), 1e5*h_fake(2), 1e5*h_fake(3), 'r', 'LineWidth', 1.3)
% % quiver3(0,0,0,1e5*h_real(1), 1e5*h_real(2), 1e5*h_real(3), 'g--', 'LineWidth', 1.3)
% quiver3(zeros(1,length(tt)),zeros(1,length(tt)),zeros(1,length(tt)), ...
%     1e5*Sun_dir_IS(1,:), 1e5*Sun_dir_IS(2,:), 1e5*Sun_dir_IS(3,:))
% title(sprintf(cspice_et2utc(t0*TU, 'C', 0)))
% subtitle('Merry XMas Enceladians!')
% legend('Saturn','FAKE','REAL','S/C')
% axis equal
% xlabel('X')
% ylabel('Y')
% zlabel('Z')