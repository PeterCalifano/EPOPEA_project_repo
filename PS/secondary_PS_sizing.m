% PS SIZING FOR DIFFERENT MISSION PHASES AND ARCHITECTURES

% Dry mass taken from ROM table without margins (because statistical analysis)

% dV_primary = [ dV_interplanetary, SOI insertion, landing ] -> bipropellant
% dV_secondary = [ SK ] -> monopropellant

%% SECONDARY SIZING
% Monopropellant 
% Small thrusters selected on the base of the mass on orbit which performs
% SK (see file word)

%% ORBITER-LANDER
clear all; close all; clc

dV = 0.275*2;
m_dry = 811 ; % total dry mass (all system does SK)
dV = sum( dV ) * 1e3 ;
Is = 235;       % MR 106E 22N thruster
m_prop = preliminary_prop_mass(dV,m_dry,Is);

%% S ORBITER + S LANDER
clear all; close all; clc

dV = 0.370*2;
m_dry =  533; % Orbiter dry mass
dV = sum( dV ) * 1e3 ;
Is = 235;       % MR 106E 22N thruster
m_prop = preliminary_prop_mass(dV,m_dry,Is);


%% NS ORBITER + S LANDER
clear all; close all; clc

dV = 0.550*2;
m_dry = 182 ; % Orbiter dry mass
dV = sum( dV ) * 1e3 ;
Is = 229;       % MR 111C 4N thruster
m_prop = preliminary_prop_mass(dV,m_dry,Is);


%% S ORBITER + MULTIPLE (2) LANDER
clear all; close all; clc

dV = 0.370*2;
m_dry = 533 ; % Orbiter dry mass
dV = sum( dV ) * 1e3 ;
Is = 235;       % MR 106E 22N thruster
m_prop = preliminary_prop_mass(dV,m_dry,Is);
