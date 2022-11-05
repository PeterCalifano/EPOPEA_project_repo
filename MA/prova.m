clear 
close 
clc

% initial & final states
RI = [1;0;0]
RF = [0;1.05;0]
VI = [0,1,0]'
VF = [-1,0,0]'

TOF = pi/2;

N = 0

M_end = 1000
n_sol = 40
n_integrator = 4
Is = 3000

[TA, time,T,m,r,flag  ] = lambLTj_plan( RI , RF , VI , VF , TOF , N  , M_end , n_sol , n_integrator,Is )

figure
plot(r.*cos(l),r.*sin(l))