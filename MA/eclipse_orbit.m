%% Illumination Conditions


cspice_kclear();

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
clear
clc
close all
mu=1.90095713928102*1e-7;
DU=238411468.296/1000; %km
TU=118760.57/(2*pi); 

mu_Sun=cspice_bodvrd('SUN','GM',1);

% Define useful constants
R_Enc = mean(cspice_bodvrd('602','RADII',3));
R_Sat = mean(cspice_bodvrd('699','RADII',3));

% Define initial time, final time, and grid of analysis:

% t0 = cspice_str2et('Jan 01 00:00:00 UTC 2051')/TU;
% tf = cspice_str2et('Jan 01 00:00:00 UTC 2055')/TU;
% tt = t0:rot:tf;

t0_try = cspice_str2et('2051 JAN 01 15:00:00.00')/TU;
%tf_try = t0_try+rot;
tf_try = cspice_str2et('2051 JAN 03 15:00:00.00')/TU;
tt_try = linspace(t0_try,tf_try,500000);

% Find the time for which Encaldus is on the x axis of Saturn Inertial
min_y = 1e10;

i_ax_Sat=(-26.73)*pi/180; 

%
i_EncSat=deg2rad(-0.009);
R_i_EncSat=[1 0 0;
        0 cos(i_EncSat) sin(i_EncSat);
        0 -sin(i_EncSat) cos(i_EncSat)];

    
for i = 1:length(tt_try)
    Sat2Enc = cspice_spkpos('602', tt_try(i)*TU, 'ECLIPJ2000', 'NONE', '699');

    %saturn orbital parameters recovery
    Sun2Sat_state=cspice_spkezr('SATURN', tt_try(i)*TU, 'ECLIPJ2000', 'NONE', 'SUN');
    r_Sat=Sun2Sat_state(1:3);
    v_Sat=Sun2Sat_state(4:6);
    
    [a_Sat,e_Sat,i_Sat,OM_Sat,om_Sat,theta_Sat] = car2kep(r_Sat,v_Sat,mu_Sun);
    %state rotation   
    %z-Rotation, om+theta
    R_an_Sat=[cos(om_Sat+theta_Sat) sin(om_Sat+theta_Sat) 0;
        -sin(om_Sat+theta_Sat) cos(om_Sat+theta_Sat) 0;
        0 0 1];
    
    %x-Rotation, i
    R_i_Sat=[1 0 0;
        0 cos(i_Sat) sin(i_Sat);
        0 -sin(i_Sat) cos(i_Sat)];
    %z-Rotation OM
    R_OM_Sat=[cos(OM_Sat) sin(OM_Sat) 0;
        -sin(OM_Sat) cos(OM_Sat) 0;
        0 0 1];
    % equator rotation
    R_eq_Sat=[1 0 0;
       0 cos(i_ax_Sat) sin(i_ax_Sat);        
       0 -sin(i_ax_Sat) cos(i_ax_Sat)];
    
    Sat2Enc=R_i_EncSat*R_an_Sat*R_i_Sat*R_OM_Sat*R_eq_Sat*Sat2Enc;
    %Sat2Enc=R_an_Sat*R_i_Sat*R_OM_Sat*Sat2Enc;
    
    if norm(Sat2Enc(2)) < min_y && Sat2Enc(1)>0
        t_start = tt_try(i);
        min_y = abs(Sat2Enc(2));
        Sat2Enc0=Sat2Enc;
    end

end

%% Eclipse computation (cit)
% Create the new time grid 
t0 = t_start;
% n_years = 3;
% n_orbits = n_years * 365 * 2; % 3 years * 365 days/year * 2 orbits/day
rot = 10000*24*3600/TU;
%tf = t0 + rot;
tf=t0+(6.034339056509405e+02*3600*24)/TU;
%tt = [linspace(t0,tf-200/TU,100),linspace(tf-199/TU,tf+2*3600/TU,20000)];
tt=linspace(tf-199/TU,tf+2*3600/TU,20000);
%sample initial state for a resonant northern L2 orbit N=4, M=11
x0_Halo=1.000062853735440;
y0_Halo=0;
z0_Halo=-0.00117884381145460;
vx0_Halo=0;
vy0_Halo=0.0168877463349484;
vz0_Halo=0;

state0_Halo=[x0_Halo,y0_Halo,z0_Halo,vx0_Halo,vy0_Halo,vz0_Halo]';

% t0=0;
% FlightDays=0.5; %days of prapagation
% tf=FlightDays*24*3600/TU; %final time of propagation
%  

%propagation - Halo
options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);
% [t_vec_Halo,state_vec_Halo]=ode113(@(t,x) CR3BP_dyn(t,x,mu),[t0,tf],state0_Halo,options_ode);

[t_vec_Halo,state_vec_Halo]=ode113(@(t,x) CR3BP_dyn(t,x,mu),tt,state0_Halo,options_ode);

state_vec_Halo=state_vec_Halo';
figure
plot3(state_vec_Halo(1,:),state_vec_Halo(2,:),state_vec_Halo(3,:));
axis equal



%initialization
x_CR3BP=state_vec_Halo(1:3,:);
x_inEnc=zeros(3,length(tt));
x_EclipSC=zeros(3,length(tt));

Sat2Enc_Fake=zeros(3,length(tt));


check_Sun = zeros(size(tt));
%Sat2Enc_array=zeros(3,length(tt));
%perc_nonvisiblity_Sun = zeros(size(tt));

for j = 1:length(tt)
    %     time_1 = tt(i-1);
    %     time_2 = tt(i);
    %     time_vec = time_1 + t_vec_Halo;
    %
    %     if tt(end) ~= time_2
    %
    %         tt(end) = time_2;
    %
    %     end
    
    % Save the time
    time_j = tt(j);
    % Compute all the relative positions between bodies
    Sat2Sun = cspice_spkpos('Sun', time_j*TU, 'ECLIPJ2000', 'NONE', '699');
    
    %saturn orbital parameters recovery
    Sun2Sat_state=cspice_spkezr('SATURN', time_j*TU, 'ECLIPJ2000', 'NONE', 'SUN');
    r_Sat=Sun2Sat_state(1:3);
    v_Sat=Sun2Sat_state(4:6);
    [a_Sat,e_Sat,i_Sat,OM_Sat,om_Sat,theta_Sat] = car2kep(r_Sat,v_Sat,mu_Sun);
    
%     Sat2Enc = cspice_spkpos('602', time_j*TU, 'ECLIPJ2000', 'NONE', '699');
%     Sat2Enc_array(:,j)=Sat2Enc;
    
    %state rotation   
    %z-Rotation, om+theta
    R_an_Sat=[cos(om_Sat+theta_Sat) sin(om_Sat+theta_Sat) 0;
        -sin(om_Sat+theta_Sat) cos(om_Sat+theta_Sat) 0;
        0 0 1]';
    %x-Rotation, i
    R_i_Sat=[1 0 0;
        0 cos(i_Sat) sin(i_Sat);
        0 -sin(i_Sat) cos(i_Sat)]';
    %z-Rotation OM
    R_OM_Sat=[cos(OM_Sat) sin(OM_Sat) 0;
        -sin(OM_Sat) cos(OM_Sat) 0;
        0 0 1]';
    % x equator - Rotation 
    R_eq_Sat=[1 0 0;
       0 cos(i_ax_Sat) sin(i_ax_Sat);        
       0 -sin(i_ax_Sat) cos(i_ax_Sat)]';
    
   % rotation to Saturn's orbital parameters 
    x_inEnc(1,j)=(x_CR3BP(1,j)+mu)*cos(tt(j)-t0) - x_CR3BP(2,j)*sin(tt(j)-t0);
    x_inEnc(2,j)=(x_CR3BP(1,j)+mu)*sin(tt(j)-t0) + x_CR3BP(2,j)*cos(tt(j)-t0);
    x_inEnc(3,j)=x_CR3BP(3,j);
    
    Sat2Enc_Fake(1,j)=1*cos(tt(j)-t0);
    Sat2Enc_Fake(2,j)=1*sin(tt(j)-t0);
    Sat2Enc_Fake(3,j)=0;
    
    Sat2Enc_Fake(:,j)=Sat2Enc_Fake(:,j)*DU;
    
    x_inEnc(:,j)=x_inEnc(:,j)*DU;
    
    x_EclipSC(:,j)=R_eq_Sat*R_OM_Sat*R_i_Sat*R_an_Sat*R_i_EncSat'*x_inEnc(:,j);
    Sat2Enc_Fake(:,j)=R_eq_Sat*R_OM_Sat*R_i_Sat*R_an_Sat*R_i_EncSat'*Sat2Enc_Fake(:,j);
    
    
    Sc2Sat = -x_EclipSC(:,j);
    
    Enc2Sun = -Sat2Enc_Fake(:,j) + Sat2Sun;
    Sc2Sun  = Sc2Sat + Sat2Sun;
    Sc2Enc  = Sc2Sat + Sat2Enc_Fake(:,j);
    
    %%% Check on Sun %%
    
    % Check if Saturn is in the way
    max_ang_Sat = atan(R_Sat/norm(Sc2Sat));
%     if max_ang_Sat > pi/4 || max_ang_Sat < 0 
%        error('Error on maximum angle of Saturn')
%     end
    ang_Sat = acos(dot(Sc2Sat,Sc2Sun)/(norm(Sc2Sat)*norm(Sc2Sun)));

    if ang_Sat < max_ang_Sat
        check_Sun(j) = check_Sun(j) + 1;
    end

    % Check if Enceladus is in the way
    max_ang_Enc = atan(R_Enc/norm(Sc2Sat));

%     if max_ang_Enc > pi/4 || max_ang_Enc < 0 
%         error('Error on maximum angle of Enceladus')
%     end

    ang_Enc = acos(dot(Sc2Enc,Sc2Sun)/(norm(Sc2Enc)*norm(Sc2Sun)));

    if ang_Enc < max_ang_Enc
        check_Sun(j) = check_Sun(j) + 2;
    end

    % Control if the vector of check_Sun is only made of zeros

%     if check_Sun == zeros(size(tt))
% 
%         % In this case the sun is always visible in the i-th orbit of 
%         % the s/c, so nothing is changed since the vector visibility_Sun 
%         % was already initialized with only ones
% 
%     else
%         
%         % The Sun is not visible
%         visibility_Sun(i) = 0;
%         
%         % Save only the instants in which the Sun is not visible
%         check_Sun_2 = check_Sun(check_Sun > 0);
%         
%         % Compute the percentage of non-visibility
%         perc_nonvisiblity_Sun = length(check_Sun_2)/length(check_Sun);
% 
%     end


end

%% Post processing

check_Sun(check_Sun>0);
tt_eclipse=zeros(1,length(check_Sun(check_Sun>0)));
tt_eclipse_str=[];
k=0;
index=zeros(1,length(check_Sun(check_Sun>0)));
for j=1:length(tt)
    if check_Sun(j)>0
        k=k+1;
        tt_eclipse(k)=(tt(j)-tt(1))*TU/(3600*24);
       index(k)=j;
    end
end


figure
plot3(x_EclipSC(1,:),x_EclipSC(2,:),x_EclipSC(3,:));
grid on
grid minor
axis equal



figure
%Saturn
scatter3(0,0,0,50);
%Enceladus
hold on
plot3(Sat2Enc_Fake(1,:),Sat2Enc_Fake(2,:),Sat2Enc_Fake(3,:),'r')
plot3(x_EclipSC(1,:),x_EclipSC(2,:),x_EclipSC(3,:),'-k')
axis equal
legend('Saturn','Enceladus','S/C')
xlabel('X')
ylabel('Y')
zlabel('Z')

figure
plot(tt,check_Sun)
xlabel('t')
ylabel('Eclipse (boolean)')
