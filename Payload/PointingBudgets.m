clear;clc;close all

%% WAC - wide remote sensing mode
h=60;
FOV_WAC=40;
Np_WAC=1024;
texp=1;
dx_WAC=100;
l_WAC=2*h*tan(deg2rad(FOV_WAC/2))*1e+3; %[m] - swath width
p_WAC=l_WAC/Np_WAC; %[m] - on ground size of a pixel

Dl=1000; %accepted drift of the centre of the image

alpha_WAC=rad2deg(atan(Dl/(h*1e+3))) %[°] pointing accuracy

p_drift=floor(dx_WAC/p_WAC); %accepted number of dislocated pixels
res=p_drift*p_WAC; %resolution
alpha_d=rad2deg(atan(res/(h*1e+3))); %drift
alpha_dot_WAC=alpha_d/texp


%% NAC - narrow remote sensing mode
FOV_NAC=0.3;
Np_NAC=1024;
texp=1;
dx_NAC=1;
[alpha_NAC,alpha_dotNAC]=pointingcam(h,FOV_NAC,Np_NAC,texp,dx_NAC)

%% TES - Wide remote sensing mode

%% TES - Narrow remote sensing mode 

%% Laser altimeter - Wide remote sensing mode

%% radar sounder - radar sounding mode

%% sample collector - passive sample collection
%implement model






