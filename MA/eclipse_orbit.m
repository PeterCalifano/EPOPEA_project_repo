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
cspice_furnsh ('spice_kernels\LATs.bsp');
cspice_furnsh ('spice_kernels\LATs.tf');
%% Propagate science orbit
clc
mu=1.90095713928102*1e-7;
DU=238411468.296/1000; %km
TU=118760.57/(2*pi); 

% Define useful constants
R_Enc = mean(cspice_bodvrd('602','RADII',3));
R_Sat = mean(cspice_bodvrd('699','RADII',3));

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

visibility_Sun = ones(size(tt));
perc_nonvisiblity_Sun = zeros(size(tt));

for i = 2:length(tt)
    time_1 = tt(i-1);
    time_2 = tt(i);    
    time_vec = time_1 + t_vec_Halo;

    if time_vec(end) ~= time_2
        time_vec(end) = time_2;
    end
    
    check_Sun = zeros(size(time_vec));

    for j = 1:length(time_vec)

        % Save the time
        time_j = time_vec(j);

        % Compute all the relative positions between bodies
        Sat2Sun = cspice_spkpos('Sun', time_j*TU, 'ECLIPJ2000', 'NONE', '699');
        Enc2Sat = cspice_spkpos('699', time_j*TU, 'ECLIPJ2000', 'NONE', '602');
        Enc2Sun = Enc2Sat + Sat2Sun;
        Sc2Enc = - state_vec_Halo(1:3,j);
        Sc2Sun  = Sc2Enc + Enc2Sun;
        Sc2Sat  = Sc2Enc + Enc2Sat;
        
        %%% Check on Sun %%%

        % Check if Saturn is in the way
        max_ang_Sat = atan(R_Sat/norm(Sc2Sat));
        if max_ang_Sat > pi/4 || max_ang_Sat < 0 
            error('Error on maximum angle of Saturn')
        end
        ang_Sat = acos(dot(Sc2Sat,Sc2Sun)/(norm(Sc2Sat)*norm(Sc2Sun)));

        if ang_Sat < max_ang_Sat
            check_Sun(j) = check_Sun(j) + 1;
        end

        % Check if Enceladus is in the way
        max_ang_Enc = atan(R_Enc/norm(Sc2Enc));

        if max_ang_Enc > pi/4 || max_ang_Enc < 0 
            error('Error on maximum angle of Enceladus')
        end

        ang_Enc = acos(dot(Sc2Enc,Sc2Sun)/(norm(Sc2Enc)*norm(Sc2Sun)));

        if ang_Enc < max_ang_Enc
            check_Sun(j) = check_Sun(j) + 2;
        end
    end

    % Control if the vector of check_Sun is only made of zeros

    if check_Sun == zeros(size(time_vec))

        % In this case the sun is always visible in the i-th orbit of 
        % the s/c, so nothing is changed since the vector visibility_Sun 
        % was already initialized with only ones

    else
        
        % The Sun is not visible
        visibility_Sun(i) = 0;
        
        % Save only the instants in which the Sun is not visible
        check_Sun_2 = check_Sun(check_Sun > 0);
        
        % Compute the percentage of non-visibility
        perc_nonvisiblity_Sun = length(check_Sun_2)/length(check_Sun);

    end


end
