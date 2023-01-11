clear;clc;close all
%DATA
g_Enc=0.113; %[m/s^2] Enceladus gravitational acceleration

g_Earth=9.81; %[m/s] earth acceleration


mb=519.2; %[kg] - lander's body mass
mlg=4*11; %[kg] - landing gear's mass

M=mb+mlg; %total lander mass
rho=2640; %[kg/m^3] - Density of the leg material (needed?)
n=4; %Number of legs

h0=2; %[m] - free fall altitude

v0=sqrt(2*g_Enc*h0); %impact velocity

ki=40e+3; %Nm - from rosetta, stiffness of a single spring
K=4*ki; %overall stiffness of the landing gear

ci=585/2; %Ns/m damper motor  
C=4*ci; %Ns/m overall damping

s=sqrt(M*v0/(K));
FcA=40; %N - amplitude of the friction force

%% simulation

t0=0;
tf=0.2; %[s] 
t_vec=linspace(t0,tf,5000);

options_ode=odeset('RelTol',1e-6,'AbsTol',1e-6);

y_vec0=[s;0];
[~,y_array]=ode113(@(t,x) singleDoFimpact(t,x,K,C,FcA,M,g_Enc),t_vec,y_vec0,options_ode);

y_array=y_array';

Fs_array=K*y_array(1,:)+C*y_array(2,:); %reaction at the ground
F_max=max(Fs_array);

figure
plot(t_vec,y_array(1,:));
hold on
plot(t_vec,s*ones(size(t_vec)),'--k')
grid minor
legend('stroke','max stroke')

figure
plot(t_vec,y_array(2,:));
grid minor
legend('velocity')

figure
plot(t_vec,Fs_array);
grid minor
legend('Force')



