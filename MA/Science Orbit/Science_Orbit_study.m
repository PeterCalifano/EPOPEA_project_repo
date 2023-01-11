%% example of science orbits
% COMMENTA STO CAZZO DI CODICE PLEASE - NO

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

% L1 determination
%Scalar potential function and gradient definition
%Symbolic manipulation
syms x y z
r1_sym= ((x+mu)^2+y^2+z^2).^0.5;
r2_sym= ((x+mu-1)^2+y^2+z^2).^0.5;

U_sym=@(x) 1/2*(x.^2+y.^2) + (1-mu)./r1_sym + mu./r2_sym+1/2*mu*(1-mu);
%dU/dx of the scalar potential function 
Ux_sym=diff(U_sym,x);

%Function handle recovery
Ux=matlabFunction(Ux_sym);

U_der_vec_x=@(x) Ux(x(1),0,0);

%zero finding
options_fzero=optimset('Display','iter','TolX',1e-15,'TolFun',1e-15);
format long

U_der_vec_x_fzero=@(x) Ux(x,0,0);
x00=1.1; %initial guess
%L2 coordinates
x_L2_fzero = fzero(U_der_vec_x_fzero,x00,options_fzero);
x_L2=[x_L2_fzero;0;0];


%CRTBP dynamics, function handle
% r1=@(t,x) ((x(1)+mu).^2+x(2).^2+x(3).^2).^0.5;
% r2=@(t,x) ((x(1)+mu-1).^2+x(2).^2+x(3).^2).^0.5;
% f=@(t,x,mu) [x(4);
%              x(5);
%              x(6);
%              2*x(5)+x(1)-(1-mu)*(x(1)+mu)/(r1(t,x)).^3 - mu/(r2(t,x)).^3*(x(1)+mu-1);
%             -2*x(4)+x(2)-(1-mu)*x(2)/r1(t,x).^3-mu/(r2(t,x)).^3*x(2);
%             -(1-mu)/(r1(t,x)).^3*x(3)-mu/(r2(t,x).^3)*x(3)];  %CRTBP dynamics
       
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
FlightDays=0.5; %days of propagation
tf=FlightDays*24*3600/TU; %final time of propagation
 

t_span = t0:(1/TU):tf;
%propagation - Halo
[t_vec_Halo,state_vec_Halo]=ode113(@(t,x) CR3BP_dyn(t,x,mu),t_span,state0_Halo,options_ode);
state_vec_Halo=state_vec_Halo';

%Position wrt Enceladus
pos_Halo_Enc=state_vec_Halo(1:3,:);
pos_Halo_Enc(1,:)=pos_Halo_Enc(1,:)+mu-1;

%in the inertial frame
%conversion of the trajectory to the inertial Enceladus Centred frame
r_vec_Halo_in=zeros(3,length(t_vec_Halo));
for k=1:length(t_vec_Halo)
    r_vec_Halo_in(1,k)=(state_vec_Halo(1,k)+mu-1)*cos(t_vec_Halo(k))-state_vec_Halo(2,k)*sin(t_vec_Halo(k));
    r_vec_Halo_in(2,k)=(state_vec_Halo(1,k)+mu-1)*sin(t_vec_Halo(k))+state_vec_Halo(2,k)*cos(t_vec_Halo(k));
    r_vec_Halo_in(3,k)=state_vec_Halo(3,k);
end

%dimension recovery
state_vec_Halo(1:3,:)=state_vec_Halo(1:3,:)*DU;
state_vec_Halo(4:6,:)=state_vec_Halo(4:6,:)*DU/TU;
r_vec_Halo_in=r_vec_Halo_in*DU; %ANTOINE - THIS IS IN THE INERTIAL FRAME 
pos_Halo_Enc=pos_Halo_Enc*DU;%ANTOINE - THIS IS IN THE ROTATING REF FRAME

%periodic resonant 3BP orbit N=9 M=35
%initial conditions
x0_P=(DU*(1-mu)-0.3628578738543753*1e+3)/DU; %km
y0_P=0.3466896463522740*1e+3/DU; %km
z0_P=0/DU; %km
vx0_P=-0.1355838999154557*1e-1/DU*TU; %km/s 
vy0_P=-0.2754409761374015*1e-1/DU*TU; %km/s
vz0_P=-0.1131541724302042/DU*TU; %km/s

state0_P=[x0_P,y0_P,z0_P,vx0_P,vy0_P,vz0_P]';
%propagation - P
tf_P=0.9896641998233017*1e+6/TU; %1 period - actually if you plot it for
%more time if gets lost (different model used in the paper)
[t_vec_P,state_vec_P]=ode113(@(t,x) CR3BP_dyn(t,x,mu),[t0 tf_P],state0_P,options_ode);
state_vec_P=state_vec_P';

% Plot Halo Orbit
Enceladus_3D(R_enc*DU,[(1-mu)*DU,0,0])
P2=plot3(state_vec_Halo(1,:),state_vec_Halo(2,:),state_vec_Halo(3,:),'k','linewidth',1.25);
P1=plot3(x_L2(1)*DU,x_L2(2)*DU,x_L2(3)*DU,'ob','markersize',5,'linewidth',1.25);
grid minor
legend([P1 P2],'$L_2$','Southern Halo orbit')

Enceladus_3D(R_enc*DU,[0,0,0])
plot3(r_vec_Halo_in(1,:),r_vec_Halo_in(2,:),r_vec_Halo_in(3,:))
grid minor 


% Plot Orbit scarabocchio 
Enceladus_3D(R_enc*DU,[(1-mu)*DU,0,0])
P1=plot3(state_vec_P(1,:)*DU,state_vec_P(2,:)*DU,state_vec_P(3,:)*DU,'b','linewidth',1);
grid minor
legend(P1,'Periodic orbit')


%% Science orbit subdivision
% propagation - half remote sensing arc (for the 3 modes)
h_RS=1; %[h] - duration of the remote sensing arc (1 of the 3 modes)
tf_RS=h_RS/2*3600/TU; 
[t_vec_RS,state_vec_RS]=ode113(@(t,x) CR3BP_dyn(t,x,mu),[t0 tf_RS],state0_Halo,options_ode);
state_vec_RS=state_vec_RS';
state_vec_RS(1:3,:)=state_vec_RS(1:3,:)*DU;
state_vec_RS(4:6,:)=state_vec_RS(4:6,:)*DU/TU;

h_CI=1; %[h] - duration of the coarse imaging arc
tf_CI=tf_RS+h_CI/2*3600/TU; 
state0_CI=[state_vec_RS(1:3,end)/DU;state_vec_RS(4:6,end)*TU/DU];
[t_vec_CI,state_vec_CI]=ode113(@(t,x) CR3BP_dyn(t,x,mu),[tf_RS tf_CI],state0_CI,options_ode);
state_vec_CI=state_vec_CI';
state_vec_CI(1:3,:)=state_vec_CI(1:3,:)*DU;
state_vec_CI(4:6,:)=state_vec_CI(4:6,:)*DU/TU;



%propagation - AOCS+sk arc/2
h_SK=3; %Number of hours dedicated to SK
tf_SK=tf_CI+h_SK/2*3600/TU; 
state0_SK=[state_vec_CI(1:3,end)/DU;state_vec_CI(4:6,end)*TU/DU];
[t_vec_SK,state_vec_SK]=ode113(@(t,x) CR3BP_dyn(t,x,mu),[tf_CI tf_SK],state0_SK,options_ode);
state_vec_SK=state_vec_SK';
state_vec_SK(1:3,:)=state_vec_SK(1:3,:)*DU;
state_vec_SK(4:6,:)=state_vec_SK(4:6,:)*DU/TU;

%propagation - communication arc
state0_else=[state_vec_SK(1:3,end)/DU;state_vec_SK(4:6,end)*TU/DU];

options_ode_event = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13,'Events',@(t,x) FirstZeroCrossing(t,x));

[t_vec_else,state_vec_else,t_e,state_e,i_e]=ode113(@(t,x) CR3BP_dyn(t,x,mu),[tf_SK tf],state0_else,options_ode_event);
state_vec_else=state_vec_else';
state_vec_else(1:3,:)=state_vec_else(1:3,:)*DU;
state_vec_else(4:6,:)=state_vec_else(4:6,:)*DU/TU;


Enceladus_3D(R_enc*DU,[(1-mu)*DU,0,0])
P1=plot3(state_vec_RS(1,:),state_vec_RS(2,:),state_vec_RS(3,:),'b','linewidth',2);
P2=plot3(state_vec_CI(1,:),state_vec_CI(2,:),state_vec_CI(3,:),'m','linewidth',2);
P3=plot3(state_vec_SK(1,:),state_vec_SK(2,:),state_vec_SK(3,:),'g','linewidth',2);
P4=plot3(state_vec_else(1,:),state_vec_else(2,:),state_vec_else(3,:),'r','linewidth',2);
plot3(state_vec_RS(1,:),-state_vec_RS(2,:),state_vec_RS(3,:),'b','linewidth',2);
plot3(state_vec_CI(1,:),-state_vec_CI(2,:),state_vec_CI(3,:),'m','linewidth',2);
plot3(state_vec_SK(1,:),-state_vec_SK(2,:),state_vec_SK(3,:),'g','linewidth',2);
plot3(state_vec_else(1,:),-state_vec_else(2,:),state_vec_else(3,:),'r','linewidth',2);
%plot3(x0_Halo*DU, y0_Halo*DU, z0_Halo*DU, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k');
grid minor
legend([P1 P2 P3 P4],'One of (1),(2),(3): 1h ','(1):1h, 0.5h per arc','SK/ADCS: 3h, 1.5h per arc','TMTC: 7h')
%% Ground Tracks 
%Ground tracks taking the state in the rotating frame
w_Enc=0; %rad/s, in the CRTBP rotating frame
[alpha, delta, lat_Halo, lon_Halo] = groundTrack(t_vec_Halo*TU, pos_Halo_Enc',90, w_Enc);

% %Ground tracks taking the state in the inertial frame
% w_Enc=1/TU; %rad/s, Enceladus centred inertial frame
% [alpha, delta, lat_Halo, lon_Halo] = groundTrack(t_vec_Halo*TU, r_vec_Halo_in',90, w_Enc);

%plot
figure
scatter(lon_Halo,lat_Halo)
xlabel('Longitude');
ylabel('Latitude');
axis equal
axis([-180 180 -90 90])
hold on
grid on
grid minor

%% Plume fly trhough
% % options_ode_flythrough = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13,'Events',@FlyThroughExit);
%  lim_lat=-65;
% % [t_vec_FT,state_vec_FT,t_lim,state_lim,i_lim]=ode113(@(t,x,lim_lat,TU,DU) CR3BP_dyn(t,x,mu),[t0 tf],state0_Halo,options_ode_flythrough,lim_lat,TU,DU);
% state_vec_FT=[];
% for k=1:length(t_vec_Halo)
%    %latitude at the k-th instant 
%    [alpha, delta, lat_t, lon_Halo] = groundTrack(t_vec_Halo(k)*TU,state_vec_Halo(:,k)',90,0);
%  if lat_t<lim_lat
%      state_vec_FT=[state_vec_FT, state_vec_Halo(:,k)];
%  else
%      break;
%  end
% 
% end
% 
% latlim=-75; %-75[deg] northest jet source latitude (according to Porco et Al and altri amici loro)
% latlim=deg2rad(latlim);
% 
% flag=0;
% for k=1:length(t_vec_Halo)
%     delta=deg2rad(lon_Halo(k)); % longitude of the s/c
%     r_enc_vec=R_enc*DU*[cos(latlim)*cos(delta);
%                         cos(latlim)*sin(delta);
%                         sin(latlim)]; %pos. vector of the northest jet source at the s/c's longitude
%     r_enc2sc=pos_Halo_Enc(:,k);
%     e_vec=r_enc2sc-r_enc_vec; %relative position s/c vs northest jet source
%     
%     r_enc_dir=r_enc_vec/norm(r_enc_vec);
%     e_dir=e_vec/norm(e_vec);
%     
%     beta=rad2deg(acos(dot(r_enc_dir,e_dir))); %angle wrt the jet
%     if beta>35 && lat_Halo(k)>latlim && flag==0 %check if the angle wrt the jet is larger than the amplitude of the jet 
%        t_exit=t_vec_Halo(k);
%        r_exit=state_vec_Halo(1:3,k);
%        r_source_exit=r_enc_vec;
%        flag=1;  
%     end
%    
% end
%     
% %propagation of the fly-through arc
% [t_vec_FT,state_vec_FT]=ode113(@(t,x) CR3BP_dyn(t,x,mu),[t0 t_exit],state0_Halo,options_ode);
% state_vec_FT=state_vec_FT';
% 
% state_vec_FT(1:3,:)=state_vec_FT(1:3,:)*DU;
% state_vec_FT(4:6,:)=state_vec_FT(4:6,:)*DU/TU;
% 
% Enceladus_3D(R_enc*DU,[(1-mu)*DU,0,0])
% P3=plot3(r_exit(1),r_exit(2),r_exit(3),'o','markersize',3,'linewidth',2)
% P2=plot3((1-mu)*DU+r_source_exit(1),r_source_exit(2),r_source_exit(3),'x','markersize',10,'linewidth',2)
% P1=plot3(state_vec_FT(1,:),state_vec_FT(2,:),state_vec_FT(3,:),'b','linewidth',2);
% %plot3(state_vec_FT(1,:),-state_vec_FT(2,:),state_vec_FT(3,:),'b','linewidth',2);
% legend([P1 P2 P3],'fly-through arc','source','exit from the plume')

flag=0;
theta_lim=-66.06; %limit latitude for plume presence (for a circular orbit with h=60km)

for k=1:length(lat_Halo)
    if lat_Halo(k)>=theta_lim && flag==0
        flag=1;
        k_lim=k;
        t_lim=t_vec_Halo(k);
    end
    
end

Time_plume=2*t_lim*TU/60 %[min] time in the plume for our science orbit

Enceladus_3D(R_enc*DU,[(1-mu)*DU,0,0])
P1=plot3(state_vec_RS(1,:),state_vec_RS(2,:),state_vec_RS(3,:),'k','linewidth',1.25);
plot3(state_vec_RS(1,:),-state_vec_RS(2,:),state_vec_RS(3,:),'k','linewidth',1.25);
P2=plot3(state_vec_Halo(1,k_lim),state_vec_Halo(2,k_lim),state_vec_Halo(3,k_lim),'or','linewidth',2,'markersize',6);
plot3(state_vec_Halo(1,k_lim),-state_vec_Halo(2,k_lim),state_vec_Halo(3,k_lim),'or','linewidth',2,'markersize',6);
grid minor
legend([P1 P2],'Halo','Fly through extremes')


%% event function
function [value,isterminal,direction] = FirstZeroCrossing(t,x) %Event function 
%to stop the integration at the intersection point with the x-z plane 

value=x(2);
isterminal=1;
direction=-1;

end

% function [value,isterminal,direction] = FlyThroughExit(t,x,lim_lat,TU,DU) %Event function 
% %to stop the integration when the s/c goes higher than a given latitude. 
% r=x(1:3)'*DU;
% time=t*TU;
% [~,~, lat_Halo,~] = groundTrack(time,r,90,0);
% value=lat_Halo-lim_lat;
% isterminal=1;
% direction=1;
% 
% end






