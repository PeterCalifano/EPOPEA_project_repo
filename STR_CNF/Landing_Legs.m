clear,clc;close all
%Script for the landing legs  study

%DATA
g_Enc=0.113; %[m/s^2] Enceladus gravitational acceleration

g_Earth=9.81; %[m/s] earth acceleration


M=519.2; %[kg] - lander mass
E= 71*10^9; %[Pa] - Young Modulus of the chosen material (Al5056)
rho=2640; %[kg/m^3] - Density of the material
n=4; %Number of legs
Ks=3; %multiplicative factor for the weight force, to account for impact and margins

%% Buckling analysis - vertical legs


%computation of the limit Pcr 
NG=3; %number of G's at impact
MaxF=M*g_Enc*Ks/n; %[N] maximum force acting on a leg, nominal operations
MaxlandF=M*NG*g_Earth/n; %[N] maximum impact force on a leg

Pcr=2*max(MaxF,MaxlandF); %[N] Margined Pcr of the single leg

%Sizing
%chosen parameters
L=1; %[m] length of the legs
th=5e-3; %[m] thickness of the leg

%solving for R, mean radius of the leg circular section
re=@(R) R+th/2;
ri=@(R) R-th/2;

I=@(R) pi*(re(R)^4-ri(R)^4)/4; %inertia moment

%Critical load equation, pinned pinned case
F_zero=@(R) Pcr-pi^2/L^2*E*I(R);

%Computation of the radius
R_guess=1e-2; %[m], initial guesss

options_fzero=optimset('Display','iter','TolX',1e-6,'TolFun',1e-6);
R_sol=fzero(F_zero,R_guess,options_fzero);

%Mass computation
M_leg=pi*(re(R_sol)^2-ri(R_sol)^2)*L*rho; %[kg], mass of 1 leg







