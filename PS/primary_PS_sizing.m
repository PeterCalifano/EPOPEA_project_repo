% PS SIZING FOR DIFFERENT MISSION PHASES AND ARCHITECTURES

% Dry mass taken from ROM table without margins (because statistical analysis)

% dV_primary = [ dV_interplanetary, SOI insertion, landing ] -> bipropellant
% dV_secondary = [ SK ] -> monopropellant

%% PRIMARY SIZING
% Bipropellant 
% 400 N Ariane Engine

%% ORBITER-LANDER
clear all; close all; clc

dV = [ 1, 3, 0.275 * 2, (0.258 + 2 * 0.007 ) ] ; % [ km/s]

m_dry = 811 ;
dV = sum( dV ) * 1e3 ;
m_prop = preliminary_prop_mass(dV,m_dry);

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

thruster.mass = 5.2 ;
thruster.number = 2 ;
thruster.chamber_pressure = 9.4e5 ;

tank.number = 3 ;
[ M_PS3, sizing_ox3, sizing_fu3, sizing_gas3 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )
tank.number = 4 ;
[ M_PS4, sizing_ox4, sizing_fu4, sizing_gas4 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )




%% S LANDER + S LANDER
clear all; close all; clc
dV = [ 1, 3, 0.370 * 2, (0.258 + 2 * 0.007 ) ] ; % [ km/s]

m_dry = 1602;  

%% NS ORBITER + S LANDER
clear all; close all; clc

m_dry = 986;


%% S ORBITER + MULTIPLE LANDER
clear all; close all; clc

m_dry = 1735;



