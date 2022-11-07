%% MISCELLANEOUS/TEST CODES - EARTH-SATURN OPTIMIZATION
clear; close all; clc

%% Load SPICE Kernels
cspice_kclear();
cspice_furnsh('..\..\spice_kernels/pck00010.tpc')
cspice_furnsh('..\..\spice_kernels/naif0012.tls')
cspice_furnsh('..\..\spice_kernels/gm_de431.tpc')
cspice_furnsh('..\..\spice_kernels/de440s.bsp')
%%
TU = 365*24*3600;
DU = 1.495978707e+8; %[Km]
date_1 = '2030-01-01 00:00:00.00 UTC';
t1 = cspice_str2et(date_1)/TU;

%% Maximum turning angle at Earth

% [max_turn_angle, hyp_params] = DetermineMaxTurnAngle(9, 'EARTH', 400)

R_Saturn = mean(cspice_bodvrd('SATURN', 'RADII', 3)); % [km]


%% JUPITER - SATURN
date_2 = '2055-01-01 00:00:00.00 UTC';
t2 = cspice_str2et(date_2)/TU;

t_v = linspace(t1,t2,100);
t_v2 = [t1:t2];
r_J = zeros(3,length(t_v));
r_S = zeros(3,length(t_v));
for i = 1:length(t_v)
    t_i = t_v(i);
    r_J(:,i) = cspice_spkpos('Jupiter Barycenter',t_i*TU,'ECLIPJ2000','NONE','SUN')/DU;
    r_S(:,i) = cspice_spkpos('Saturn Barycenter',t_i*TU,'ECLIPJ2000','NONE','SUN')/DU;
end

figure
scatter3(0,0,0,50,'yellow','filled')
hold on
grid on
xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')
axis equal
plot3(r_J(1,:),r_J(2,:),r_J(3,:),'k','linewidth',1.5)
plot3(r_S(1,:),r_S(2,:),r_S(3,:),'k','linewidth',1.5)


r_J2 = zeros(3,length(t_v2));
r_S2 = zeros(3,length(t_v2));
for i = 1:length(t_v2)
    t_i = t_v2(i);
    r_J2(:,i) = cspice_spkpos('Jupiter Barycenter',t_i*TU,'ECLIPJ2000','NONE','SUN')/DU;
    r_S2(:,i) = cspice_spkpos('Saturn Barycenter',t_i*TU,'ECLIPJ2000','NONE','SUN')/DU;
end

j = scatter3(1,1,1,20);
s = scatter3(1,1,1,20);
for i = 1:length(t_v2)
    delete(j)
    delete(s)
    j = scatter3(r_J2(1,i),r_J2(2,i),r_J2(3,i),30,'r','filled');
    s = scatter3(r_S2(1,i),r_S2(2,i),r_S2(3,i),30,'blue','filled');
    year = 2030 + i -1;
    title(['Year = ',num2str(year)])
    drawnow
    pause(1)
end

%% VENUS - EARTH
clc; close all; clear;
TU = 365*24*3600;
DU = 1.495978707e+8; %[Km]
date_1 = '2030-01-01 00:00:00.00 UTC';
t1 = cspice_str2et(date_1)/TU;
t1 = t1*365;
t2 = t1 + 584;

t_v = linspace(t1,t2,1000);
t_v2 = [t1:30:t2];
r_V = zeros(3,length(t_v));
r_E = zeros(3,length(t_v));
for i = 1:length(t_v)
    t_i = t_v(i);
    r_V(:,i) = cspice_spkpos('Venus',t_i*TU/365,'ECLIPJ2000','NONE','SUN')/DU;
    r_E(:,i) = cspice_spkpos('Earth',t_i*TU/365,'ECLIPJ2000','NONE','SUN')/DU;
end

figure
scatter3(0,0,0,50,'yellow','filled')
hold on
grid on
xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')
axis equal
plot3(r_V(1,:),r_V(2,:),r_V(3,:),'k','linewidth',1.5)
plot3(r_E(1,:),r_E(2,:),r_E(3,:),'k','linewidth',1.5)

r_V2 = zeros(3,length(t_v2));
r_E2 = zeros(3,length(t_v2));
for i = 1:length(t_v2)
    t_i = t_v2(i);
    r_V2(:,i) = cspice_spkpos('Venus',t_i*TU/365,'ECLIPJ2000','NONE','SUN')/DU;
    r_E2(:,i) = cspice_spkpos('Earth',t_i*TU/365,'ECLIPJ2000','NONE','SUN')/DU;
end

j = scatter3(1,1,1,20);
s = scatter3(1,1,1,20);
for i = 1:length(t_v2)
    delete(j)
    delete(s)
    j = scatter3(r_V2(1,i),r_V2(2,i),r_V2(3,i),30,'r','filled');
    s = scatter3(r_E2(1,i),r_E2(2,i),r_E2(3,i),30,'blue','filled');
    if i <= 12
        year = 2030;
        month = i;
    elseif i > 12
        year = 2031;
        month = i - 12;
    end
    title(['Year = ',num2str(year),'  Month = ',num2str(month)])
    view([0,0,1])
    drawnow
    pause()
end

%% OPTIMIZATION
R_S = cspice_bodvrd('Saturn','RADII',3);
R_S = R_S(1);
N = 3;
Ra_target = 200*R_S;
Rp_target = 20*R_S;
planets = {'Earth','Venus','Earth','Jupiter Barycenter','Saturn Barycenter'};
t1 = cspice_str2et('2030-01-01 00:00:00.00 UTC');
TU = 365*24*3600; % normalization in years
t1 = t1/TU;
guess_0 = [t1,2/12,2/12,4,5];
%%
clc
A = [];
b = [];
DV = fmincon(@(var) objfun_EarthSaturntransfer(var,N,planets,Ra_target,Rp_target,TU),guess_0,A,b);



