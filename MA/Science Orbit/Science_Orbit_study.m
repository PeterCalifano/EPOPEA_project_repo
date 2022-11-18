%% example of a science orbit plot

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(0,'defaultAxesFontSize', 15)
clear;clc;close all

%units
mu=1.90095713928102*1e-7;
DU=238411468.296/1000; %km
TU=118760.57; %s

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

%sample initial state for a resonant southern L2 orbit N=4, M=11
x0=1.000062853735440;
y0=0;
z0=-0.00117884381145460;
vx0=0;
vy0=0.0168877463349484;
vz0=0;

state0=[x0,y0,z0,vx0,vy0,vz0]';

t0=0;
FlightDays=30; %days of prapagation
tf=FlightDays*24*3600/TU; %final time of propagation

%propagation
[t_vec,state_vec]=ode113(@(t,x) f(t,x,mu),[t0 tf],state0,options_ode);
state_vec=state_vec';

%plot
Enceladus_3D_Adim(R_enc,[1-mu,0,0])
plot3(state_vec(1,:),state_vec(2,:),state_vec(3,:),'k')




