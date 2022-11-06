% PS SIZING FOR DIFFERENT MISSION PHASES AND ARCHITECTURES

% Dry mass taken from ROM table without margins (because statistical analysis)

% dV_primary = [ dV_interplanetary, SOI insertion, landing ] -> bipropellant
% dV_secondary = [ SK ] -> monopropellant

%% PRIMARY SIZING
% Bipropellant 
% 400 N Ariane Engine

% Pressurant storage pressure taken from tanks datasheets
% https://www.northropgrumman.com/space/pressurant-tanks-data-sheets-sorted-by-volume/

% Tanks material data preliminarly taken from slides


%% ORBITER-LANDER + SM
clear all; close all; clc

% dV = [interplanetary, capture, Moon tour, landing] + stochastic margins;
dV = [ 0.5*2, 0.9, 1.1*2, 0.25+0.007*2 ] ; % [ km/s]

m_dry = 811 ;
dV = sum( dV ) * 1e3 ;
Is = 321;       % Ariane engine
m_prop = preliminary_prop_mass(dV,m_dry,Is);
m_wet = m_prop+m_dry; % PRIMARY wet

% % Bipropellant
% oxidizer.rho = 1443 ;
% oxidizer.OFratio = 1.65 ;
% oxidizer.prop_mass = m_prop ;
% 
% fuel.rho = 900 ;
% 
% pressurant.feed_system = 'PressureRegulated' ;
% pressurant.pressure = 300e5 ;
% pressurant.gamma = 1.667 ;
% pressurant.R = 2077.3 ;
% 
% % SPHERICAL CASE
% tank.type = 'Spherical' ;
% tank.rho = 2780;
% tank.sigma = 950e6 ;
% 
% thruster.mass = 4.30 ;
% thruster.number = 2 ;
% thruster.chamber_pressure = 10.35e5 ;
% 
% tank.number = 3 ;       % 1 for pressurant, 2 for ox/fuel
% [ M_PS3, sizing_ox3, sizing_fu3, sizing_gas3 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )
% tank.number = 4 ;       % 2 for pressurant, 2 for ox/fuel
% [ M_PS4, sizing_ox4, sizing_fu4 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )
% 
% % CYLINDRICAL CASE
% % h = 
% tank.type = 'Cylindrical' ;
% 

%% S ORBITER + S LANDER
clear all; close all; clc

m_dry = 1279; % total dry mass (orbiter+lander) 
Is_int = 321;       % Ariane engine
% INTERPLANETARY+INSERTION+MOON
% dV = [interplanetary, capture, Moon tour] + stochastic margins;
dV_int = [ 0.5*2, 0.9, 1.1*2 ] ; % [ km/s]
 
dV_int = sum( dV_int ) * 1e3 ;
m_prop_int = preliminary_prop_mass(dV_int,m_dry,Is_int);

% LANDING
m_dry_landing = 746 ; % lander dry mass
Is_landing = 228;
dV_landing = 0.25+0.007*2;
% dV = landing + stochastic margins
dV_landing = dV_landing*1e+3;
m_prop_landing = preliminary_prop_mass(dV_landing,m_dry_landing,Is_landing);

m_wet = m_dry + m_prop_int + m_prop_landing;
% % Bipropellant
% oxidizer.rho = 1443 ;
% oxidizer.OFratio = 1.65 ;
% oxidizer.prop_mass = m_prop ;
% 
% fuel.rho = 900 ;
% 
% pressurant.feed_system = 'PressureRegulated' ;
% pressurant.pressure = 300e5 ;
% pressurant.gamma = 1.667 ;
% pressurant.R = 2077.3 ;
% 
% % SPHERICAL CASE
% tank.type = 'Spherical' ;
% tank.rho = 2780;
% tank.sigma = 950e6 ;
% 
% thruster.mass = 4.30 ;
% thruster.number = 2 ;
% thruster.chamber_pressure = 10.35e5 ;
% 
% tank.number = 3 ;       % 1 for pressurant, 2 for ox/fuel
% [ M_PS3, sizing_ox3, sizing_fu3, sizing_gas3 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )
% tank.number = 4 ;       % 2 for pressurant, 2 for ox/fuel
% [ M_PS4, sizing_ox4, sizing_fu4 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )

%% NS ORBITER + S LANDER
clear all; close all; clc

m_dry = 986;        % total dry (orbiter + lander) 
Is_int = 321;       % Ariane engine
% INTERPLANETARY+INSERTION+MOON
% dV = [interplanetary, capture, Moon tour] + stochastic margins;
dV_int = [ 0.5*2, 0.9, 1.1*2 ] ; % [ km/s]
 
dV_int = sum( dV_int ) * 1e3 ;
m_prop_int = preliminary_prop_mass(dV_int,m_dry,Is_int);

% LANDING
m_dry_landing = 746 ; % lander dry mass
Is_landing = 228;
dV_landing = 0.25+0.007*2;
% dV = landing + stochastic margins
dV_landing = dV_landing*1e+3;
m_prop_landing = preliminary_prop_mass(dV_landing,m_dry_landing,Is_landing);

m_wet = m_dry + m_prop_int + m_prop_landing;
% % Bipropellant
% oxidizer.rho = 1443 ;
% oxidizer.OFratio = 1.65 ;
% oxidizer.prop_mass = m_prop ;
% 
% fuel.rho = 900 ;
% 
% pressurant.feed_system = 'PressureRegulated' ;
% pressurant.pressure = 300e5 ;
% pressurant.gamma = 1.667 ;
% pressurant.R = 2077.3 ;
% 
% % SPHERICAL CASE
% tank.type = 'Spherical' ;
% tank.rho = 2780;
% tank.sigma = 950e6 ;
% 
% thruster.mass = 4.30 ;
% thruster.number = 2 ;
% thruster.chamber_pressure = 10.35e5 ;
% 
% tank.number = 3 ;       % 1 for pressurant, 2 for ox/fuel
% [ M_PS3, sizing_ox3, sizing_fu3, sizing_gas3 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )
% tank.number = 4 ;       % 2 for pressurant, 2 for ox/fuel
% [ M_PS4, sizing_ox4, sizing_fu4 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )
% 

%% S ORBITER + MULTIPLE (2) LANDER
clear all; close all; clc

m_dry = 1867;       % orbiter + 2 landers
Is_int = 321;       % Ariane engine
% INTERPLANETARY+INSERTION+MOON
% dV = [interplanetary, capture, Moon tour] + stochastic margins;
dV_int = [ 0.5*2, 0.9, 1.1*2 ] ; % [ km/s]
 
dV_int = sum( dV_int ) * 1e3 ;
m_prop_int = preliminary_prop_mass(dV_int,m_dry,Is_int);

% LANDING
m_dry_landing = 1594 ; % landers dry mass (considers both landers)
Is_landing = 228;
dV_landing = 0.25+0.007*2;
% dV = landing + stochastic margins
dV_landing = 2*dV_landing*1e+3;
m_prop_landing = preliminary_prop_mass(dV_landing,m_dry_landing,Is_landing);

m_wet = m_dry + m_prop_int + m_prop_landing;
% 
% % Bipropellant
% oxidizer.rho = 1443 ;
% oxidizer.OFratio = 1.65 ;
% oxidizer.prop_mass = m_prop ;
% 
% fuel.rho = 900 ;
% 
% pressurant.feed_system = 'PressureRegulated' ;
% pressurant.pressure = 300e5 ;
% pressurant.gamma = 1.667 ;
% pressurant.R = 2077.3 ;
% 
% % SPHERICAL CASE
% tank.type = 'Spherical' ;
% tank.rho = 2780;
% tank.sigma = 950e6 ;
% 
% thruster.mass = 4.30 ;
% thruster.number = 2 ;
% thruster.chamber_pressure = 10.35e5 ;
% 
% tank.number = 3 ;       % 1 for pressurant, 2 for ox/fuel
% [ M_PS3, sizing_ox3, sizing_fu3, sizing_gas3 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )
% tank.number = 4 ;       % 2 for pressurant, 2 for ox/fuel
% [ M_PS4, sizing_ox4, sizing_fu4 ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )




