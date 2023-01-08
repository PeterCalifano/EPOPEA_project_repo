%% Propagation of the whole Halo
clear,clc;close all

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(0,'defaultAxesFontSize', 13)
clear;clc;close all

%units
mu=1.90095713928102*1e-7;
DU=238411468.296/1000; %km
TU=118760.57/(2*pi); 

R_enc=252.1; %km, Enceladus mean radius
R_enc=R_enc/DU;

options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);

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
 

t_span = t0:(1/TU):tf;
%propagation - Halo
[t_vec_Halo,state_vec_Halo]=ode113(@(t,x) CR3BP_dyn(t,x,mu),t_span,state0_Halo,options_ode);
state_vec_Halo=state_vec_Halo';

%Position wrt Enceladus
pos_Halo_Enc=state_vec_Halo(1:3,:);
pos_Halo_Enc(1,:)=pos_Halo_Enc(1,:)+mu-1;

%dimension recovery
state_vec_Halo(1:3,:)=state_vec_Halo(1:3,:)*DU;
state_vec_Halo(4:6,:)=state_vec_Halo(4:6,:)*DU/TU;
pos_Halo_Enc=pos_Halo_Enc*DU;


% science arc recovery
h_RS=2; %[h] - duration of the whole science arc
tf_RS=h_RS/2*3600/TU; 
t_span_RS = t0:(1/TU):tf_RS;


r_RS_Enc=pos_Halo_Enc(:,1:length(t_span_RS));
r_RS_Enc_sym=[r_RS_Enc(1,end:-1:2);-r_RS_Enc(2,end:-1:2);r_RS_Enc(3,end:-1:2)];
r_RS_Enc=[r_RS_Enc_sym,r_RS_Enc];


t_vec_RS=[-t_span_RS(end:-1:2),t_span_RS];

%plot of the trajectories, around Enceladus 
Enceladus_3D(R_enc*DU,[0,0,0])
plot3(pos_Halo_Enc(1,:),pos_Halo_Enc(2,:),pos_Halo_Enc(3,:),'k');
plot3(r_RS_Enc(1,:),r_RS_Enc(2,:),r_RS_Enc(3,:),'r');

%% Coarse Imaging mode
Max_Pics=1000000000;
w_Enc=0; %rad/s, in the CRTBP rotating frame
FoV_WAC=40; %deg
data_WAC=0.0042; %Gb
Overlap_WAC=50; % % of overlap  

[track_length_CI,Nptot_CI,Data_tot_CI,rev_CI,Nporb_CI,Data_Orb_CI,left_swath_CI,right_swath_CI,print] = GroundImaging(t_vec_RS*TU,r_RS_Enc,R_enc*DU,270,w_Enc,FoV_WAC,data_WAC,Overlap_WAC,Max_Pics);






%% Ground Track

%Ground tracks taking the state in the rotating frame

[alpha, delta, lat_Halo, lon_Halo] = groundTrack(t_vec_Halo*TU, pos_Halo_Enc',270, w_Enc);

[alpha_RS, delta_RS, lat_RS, lon_RS] = groundTrack(t_vec_RS*TU, r_RS_Enc',270, w_Enc);

%plot
figure
P1=scatter(lon_Halo,lat_Halo,'.')
hold on
P2=scatter(lon_RS,lat_RS,'.')
P3=plot(left_swath_CI(1,:),left_swath_CI(2,:),'k','linewidth',1.25)
plot(right_swath_CI(1,:),right_swath_CI(2,:),'k','linewidth',1.25)
plot(print(1,:),print(2,:),'linewidth',3); 
xlabel('Longitude');
ylabel('Latitude');
axis equal
axis([-180 180 -90 90])
hold on
grid on
grid minor
legend([P1 P2 P3],'Halo','Remote sensing arc','CI swath')


%next steps:
% - compute track length on the ground
% - compute the number of pictures to cover it
% - start thinking about swaths


