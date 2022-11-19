%% example of science orbits

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(0,'defaultAxesFontSize', 15)
clear;clc;close all

%units
mu=1.90095713928102*1e-7;
DU=238411468.296/1000; %km
TU=118760.57/(2*pi); 

R_enc=252.1; %km, Enceladus mean radius



R_enc=R_enc/DU;

%CRTBP dynamics, function handle
r1=@(t,x) ((x(1)+mu).^2+x(2).^2+x(3).^2).^0.5;
r2=@(t,x) ((x(1)+mu-1).^2+x(2).^2+x(3).^2).^0.5;
f=@(t,x,mu) [x(4);
             x(5);
             x(6);
             2*x(5)+x(1)-(1-mu)*(x(1)+mu)/(r1(t,x)).^3 - mu/(r2(t,x)).^3*(x(1)+mu-1);
            -2*x(4)+x(2)-(1-mu)*x(2)/r1(t,x).^3-mu/(r2(t,x)).^3*x(2);
            -(1-mu)/(r1(t,x)).^3*x(3)-mu/(r2(t,x).^3)*x(3)];  %CRTBP dynamics
       
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
%FlightDays=2; %days of prapagation
%tf=FlightDays*24*3600/TU; %final time of propagation
tf=1*3600/TU;  


%propagation - Halo
[t_vec_Halo,state_vec_Halo]=ode113(@(t,x) f(t,x,mu),[t0 tf],state0_Halo,options_ode);
state_vec_Halo=state_vec_Halo';
state_vec_Halo(1:3,:)=state_vec_Halo(1:3,:)*DU;
state_vec_Halo(4:6,:)=state_vec_Halo(4:6,:)*DU/TU;

%plot
% Enceladus_3D_Adim(R_enc,[1-mu,0,0])
% plot3(state_vec_Halo(1,:),state_vec_Halo(2,:),state_vec_Halo(3,:),'k')




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
[t_vec_P,state_vec_P]=ode113(@(t,x) f(t,x,mu),[t0 tf_P],state0_P,options_ode);
state_vec_P=state_vec_P';



Enceladus_3D(R_enc*DU,[(1-mu)*DU,0,0])
plot3(state_vec_Halo(1,:),state_vec_Halo(2,:),state_vec_Halo(3,:),'k')
%plot3(state_vec_P(1,:)*DU,state_vec_P(2,:)*DU,state_vec_P(3,:)*DU,'b')

% %% ground tracks
% w_Enc=1/TU;
% 
% 
% %conversion of the trajectory to the inertial Enceladus Centred frame
% for k=1:length(t_vec_guess)
%     x_Halo_GT(k)=((1,k)+mu)*cos(t_vec_guess(k))-xx(2,k)*sin(t_vec_guess(k));
%     y_Halo_GT(k)=(xx(1,k)+mu)*sin(t_vec_guess(k))+xx(2,k)*cos(t_vec_guess(k));
% end
% om_em=2.66186135e-6; %E-M angular velocity
% figure
% plot(X1,Y1,'r','linewidth',1.5);
% hold on
% plot(0,0,'xb','Markersize',4,'linewidth',3)
% t_plot=t_vec_guess*TU*24*3600;
% plot(cos(om_em*t_plot),sin(om_em*t_plot),'--k','linewidth',1.5);
% grid on
% grid minor
% axis equal
% xlabel('x')
% ylabel('y')
% legend('Guess trajectory','Earth','Motion of the Moon','location','Northwest')
% title('Guess trajectory, non dimensional inertial frame')
% 
% 
% 
% 
% 
% [alpha, delta, lat_Halo, lon_Halo] = groundTrack(t_vec_Halo*TU, state_vec_Halo_GT,90, w_Enc);
% 
% figure
% plot(lon_Halo,lat_Halo,'k')
% xlabel('Longitude');
% ylabel('Latitude');
% hold on
% grid on
% grid minor





