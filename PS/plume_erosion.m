clear all, close all, clc

set( 0, 'defaulttextinterpreter', 'latex' )
set( 0, 'defaultlegendinterpreter', 'latex' )
set( groot,'defaultaxesfontSize', 17 )

% ESTIMATE THRUSTER PLUME EROSION AND CONTAMINATION
% MR107T 100 N thruster (hydrazine)

%% Mass rate deposition

m_lander = 943.4437 ; % [ kg ] - Lander total mass (from architecture_sizing.m)
g = 0.1 ; % [ m/s^2 ] - Enceladus gravitational acceleration
T = m_lander * g ; % [ N ] - Thrust needed during landing (assumed to be equal to lander weight)

Isp = 228 ; % [ s ] - Specific impulse of MR107T thruster
g0 = 9.81 ; % [ m/s ] - Gravitational acceleration on Earth at sea level
mdot = T / ( Isp * g0 ) ; % [ kg/s ] - Mass flow raterequired to generate thrust T

alpha = 45 ; % [ deg ] - Half-cone angle representing exhaust plume expansion -> height and ardius of the cone are equal
h_vect = linspace( 2, 50, 50 ) ; % [ m ] - Height from ground at which we turn off the thrusters
v = 0.75 ; % [ m / s ] - Velocity at landing (assumed from literature)

m_fluence = zeros( length( h_vect ), 1 ) ; % Initialize vector

for k = 1 : length( h_vect )

    h = h_vect( k ) ;
    m_fluence( k ) = T / ( Isp * g0 * pi * h * v ) ; % [ kg / m^2 ] - Derivative of mass flux within the cone over time

end

m_fluence = m_fluence * 1e3 ; % [ g/m^2]

scale = 0.1 ; % Scaling factor (account for percentage of ammonia in thruster plume) - from literature (Phoenix lander)

% Plot the mass fluence as function of cutoff height
figure()
plot( h_vect, m_fluence, 'Linewidth', 2 )
hold, grid on
plot( h_vect, m_fluence * scale, 'Linewidth', 2 )
xlabel( 'Cutoff height [m]' ), ylabel( 'Mass fluence [$g/m^2$]' )
title( 'Exhaust plume deposition' )
legend( 'Exhaust plume deposition', 'Ammonia deposition', 'location', 'best' )

%% Plume erosion



