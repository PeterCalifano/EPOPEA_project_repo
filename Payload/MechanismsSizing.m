clear;clc;close all



% Sample handling mechanism sizing
rho=2100; %[kg/m^3], aluminum's density 

%radii of the ring
ri=1; %[m] from conf
re=1.01;%[m] to be assumed

%thicknesses of the ring
h_PL_mounting=0.01; %[m] from conf
h_box=0.01; %[m], to be sized, for now assumed
h=h_box+h_PL_mounting;

b=0.5; %[Nms/rad] - linear damping coefficient, totally invented
theta_dot=2*pi/(3600); %[rad/s] to do a full lap in 1 hour
theta_ddot=0.2; %[rad/s^2]

eta_mot=0.7; %[invented]

m=pi*(re^2-ri^2)*h*rho; %[kg]

J=1/2*m*(ri^2+re^2); %[kg m^2], empty cilinder


%sizing torques
T_motion=b*theta_dot; %[Nm]
T_startup=J*theta_ddot; %[Nm]  
T=max(T_startup,T_motion); 



SHM.MassB=2*m; %Kg, 100%margin
SHM.PowerB=2*T*theta_dot/eta_mot %[W], 100% margin


