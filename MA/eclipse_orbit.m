%% Illumination Conditions
% SATURN USED 

clear
clc
close all

cspice_kclear()

% CHANGE WITH YOUR SCRIPT DIR: IT MUST CONTAIN ALSO THE KERNEL FOLDER
% cd('C:\Users\26ele\Desktop\POLI 2\ASMAD\PROJECT\Matlab Codes')

% kernelpool = fullfile('EPOPEA_metakernel.tm'); % for Windows
% cspice_furnsh(kernelpool);

% Define Kernel
cspice_furnsh('spice_kernels/naif0012.tls')
cspice_furnsh('spice_kernels/de440s.bsp')
cspice_furnsh('spice_kernels/sat441.bsp')
cspice_furnsh('spice_kernels/pck00010_mod.tpc')
cspice_furnsh('spice_kernels/gm_de431.tpc')
cspice_furnsh('spice_kernels/dss.bsp')
cspice_furnsh('spice_kernels/dss.tf')
cspice_furnsh('spice_kernels/pck00010.tpc')
cspice_furnsh('spice_kernels/plu058.bsp')
%% Propagate science orbit
clc
mu=1.90095713928102*1e-7;
DU=238411468.296/1000; %km
TU=118760.57/(2*pi); 

% Define initial time, final time, and grid of analysis:
rot = 0.5*24*3600/TU;
t0 = cspice_str2et('Jan 01 00:00:00 UTC 2051')/TU;
tf = cspice_str2et('Jan 01 00:00:00 UTC 2055')/TU;
tt = t0:rot:tf;


%sample initial state for a resonant northern L2 orbit N=4, M=11
x0_Halo=1.000062853735440;
y0_Halo=0;
z0_Halo=-0.00117884381145460;
vx0_Halo=0;
vy0_Halo=0.0168877463349484;
vz0_Halo=0;

state0_Halo=[x0_Halo,y0_Halo,z0_Halo,vx0_Halo,vy0_Halo,vz0_Halo]';

t0=0;
FlightDays=0.5; %days of prapagation
tf=FlightDays*24*3600/TU; %final time of propagation
 
%propagation - Halo
options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);
[t_vec_Halo,state_vec_Halo]=ode113(@(t,x) CR3BP_dyn(t,x,mu),[t0,tf],state0_Halo,options_ode);
state_vec_Halo=state_vec_Halo';
state_vec_Halo(1:3,:)=state_vec_Halo(1:3,:)*DU;
state_vec_Halo(4:6,:)=state_vec_Halo(4:6,:)*DU/TU;

for i = 2:length(tt)
    time_1 = tt(i-1);
    time_2 = tt(i);    
    time_vec = time_1 + t_vec_Halo;

    if time_vec(end) ~= time_2
        time_vec(end) = time_2;
    end

    for j = 1:length(time_vec)
        time_j = time_vec(j);
        Sat2Sun = cspice_spkpos('Sun', time_j*TU, 'ECLIPJ2000', 'NONE', '699');
        Enc2Sat = cspice_spkpos('699', time_j*TU, 'ECLIPJ2000', 'NONE', '602');
        Enc2Sun = Sat2Sun + Enc2Sat;
        Sc2Enc = state_vec_Halo(1:3,j);
        Sc2Sun  = Sc2Enc + Enc2Sun;

    end


end
%%




ett_obs = tt*TU;




max_el90 = zeros(1, length(ett_obs));
max_el80 = zeros(1, length(ett_obs));
max_el60 = zeros(1, length(ett_obs));
max_el0 = zeros(1, length(ett_obs));
for i=1:length(ett_obs)
    if i == length(ett_obs)
        continue
    end
    t1 = ett_obs(i);
    t2 = ett_obs(i+1);
    tspan = t1:60:t2;

    state90 = cspice_spkpos('SUN', tspan, 'LAT90_TOPO', 'NONE', 'LAT90');
    state80 = cspice_spkpos('SUN', tspan, 'LAT80_TOPO', 'NONE', 'LAT80');
    state60 = cspice_spkpos('SUN', tspan, 'LAT60_TOPO', 'NONE', 'LAT60');
    state0 = cspice_spkpos('SUN', tspan, 'LAT0_TOPO', 'NONE', 'LAT0');
    [rho90, az90, el90] = cspice_reclat(state90);
    [rho80, az80, el80] = cspice_reclat(state80);
    [rho60, az60, el60] = cspice_reclat(state60);
    [rho0, az0, el0] = cspice_reclat(state0);
    max_el90(i) = max(rad2deg(el90));
    max_el80(i) = max(rad2deg(el80));
    max_el60(i) = max(rad2deg(el60));
    max_el0(i) = max(rad2deg(el0));
end

%% Plots
close all; clc

figure;
time = ett_obs/86400;
shift = date2jd([2000,01,01,12,00,00])-date2jd([0000,01,01,12,00,00]);
time = time + shift;
plot(time, max_el90,'color', '#0072BD', 'LineWidth', 2);
hold on
plot(time, max_el80,'color', '#D95319', 'LineWidth', 2);
plot(time, max_el60,'color', '#EDB120', 'LineWidth', 2);
plot(time, max_el0,'color', '#77AC30', 'LineWidth', 2);
grid on
grid minor
ax = gca;
ax.XTick = time(1):2*365:time(end);
datetick('x', 'yyyy mm dd', 'keeplimits', 'keepticks');
xtickangle(45);
xlim([time(1) time(end)]);
ylim([0 90]);

ylabel('Maximum Elevation [deg]')
s = legend('0° long, -90° lat', '0° long, -80° lat', '0° long, -60° lat', '0° long, 0° lat');
s.Location = 'best';
