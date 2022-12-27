clear,clc;close all
%Script for the landing legs  study

%DATA
g_Enc=0.113; %[m/s^2] Enceladus gravitational acceleration

M=462.4; %[kg] - lander mass
E= 71*10^9; %[Pa] - Young Modulus of the chosen material (Al5056)
rho=2640; %[kg/m^3] - Density of the material
n=4; %Number of legs
Ks=3; %multiplicative factor for the weight force, to account for impact and margins

%% Buckling analysis - vertical legs


%computation of the limit Pcr 
MaxF=M*g_Enc*Ks/n; %maximum force acting on a leg
Pcr=1.5*MaxF; %Margined pcr of the single leg

%Sizing
%chosen parameters
L=1; %[m] length of the legs
th=2e-3; %[m] thickness of the leg

%solving for R, mean radius of the leg circular section
re=@(R) R+th/2;
ri=@(R) R-th/2;

I=@(R) pi*(re(R)^4-ri(R)^4)/4;

%Critical load equation, pinned pinned case
F_zero=@(R) Pcr-pi^2/L^2*E*I(R);

%Computation of the radius
R_guess=1e-2; %[m], initial guesss

options_fzero=optimset('Display','iter','TolX',1e-6,'TolFun',1e-6);
R_sol=fzero(F_zero,R_guess,options_fzero);

%Mass computation
M_leg=pi*(re(R_sol)^2-ri(R_sol)^2)*L*rho; %[kg], mass of 1 leg







