% PS SIZING FOR DIFFERENT MISSION PHASES AND ARCHITECTURES

% Dry mass taken from ROM table without margins (because statistical analysis)

% dV_primary = [ dV_interplanetary, SOI insertion, landing ] -> bipropellant
% dV_secondary = [ SK ] -> monopropellant


%% ORBITER-LANDER
clear all
close all
clc

dV = [ 1, 3, 0.275 * 2, (0.258 + 2 * 0.007 ) ] ; % [ km/s]

m_dry = 881 ;
dV = sum( dV ) * 1e3 ;
m_prop = preliminary_prop_mass(dV,m_dry)

% Bipropellant
oxidizer.rho = 1443 ;
oxidizer.OFratio = 0.85 ;
oxidizer.prop_mass = m_prop ;

fuel.rho = 900 ;

pressurant.feed_system = 'PressureRegulated' ;
pressurant.pressure = 300e5 ;
pressurant.gamma = 1.667 ;
pressurant.R = 2077.3 ;

tank.type = 'Spherical' ;
tank.rho = 2780 ;
tank.sigma = 950e6 ;
tank.number = 3 ;

thruster.mass = 5.2 ;
thruster.number = 2 ;
thruster.chamber_pressure = 9.4e5 ;

[ M_PS, sizing_ox, sizing_fu, sizing_gas ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )

%% S LANDER + S LANDER
clear all
close all
clc
dV = [ 1, 3, 0.370 * 2, (0.258 + 2 * 0.007 ) ] ; % [ km/s]

m_dry = 1309 ;
dV = sum( dV ) * 1e3 ;
m_prop = preliminary_prop_mass(dV,m_dry)