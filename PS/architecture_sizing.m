% Architectures sizing considering a specific tank configuration

%% SAMPLING ORBITER - SAMPLING LANDER
close all, clear all, clc

%%%%%%%%
% GENERAL DATA
pressurant.R = 2077.3 ;
pressurant.gamma = 1.667 ;

tank.rho = 2780 ;
tank.sigma = 950e6 ;
%%%%%%%%

% BIPROPELLANT - 400N BI-PROPELLANT APOGEE THRUSTER
dV_int = 30 + 1.05*0.5e+3 + 1.05*0.9e+3 + 2*1.1e+3 + 10;  % plus MAR-DV-010 / MAR-DV-020 / MAR-DV-080 / MAR-DV-100
Isp_int = 321; m_dry_int = 1279; % Orbiter + lander
m_prop_int = preliminary_prop_mass(dV_int,m_dry_int,Isp_int);

%m_prop_int = 3491.9;                    % from preliminary sizing

oxidizer.rho = 1443 ; 
oxidizer.OFratio = 1.65 ;
oxidizer.prop_mass = m_prop_int ;
oxidizer.type = 'Cylindrical' ;

fuel.rho = 900 ;

pressurant.feed_system = 'PressureRegulated' ;
pressurant.pressure = 300e5 ; 
pressurant.type = 'Spherical' ;

tank.ratio = 1 ;
tank.number = 3 ;

thruster.mass = 4.30 ;
thruster.number = 2 ;
thruster.chamber_pressure = 10.35e5 ;

[ M_PS_SOSL_INT, sizing_ox_SOSL_INT, sizing_fu_SOSL_INT, sizing_gas_SOSL_INT ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )

% MONOPROPELLANT - LANDING - MR107T 110 N
dV_land = 258*2;                            % Additional MAR-DV-020
Isp_land = 228; m_dry_land = 746;           % Sampling lander
m_prop_land = preliminary_prop_mass(dV_land,m_dry_land,Isp_land);

%m_prop_land = 95.33;

propellant.mass = m_prop_land ;
propellant.rho = 1e3 ;
propellant.type = 'Cylindrical' ;

pressurant.feed_system = 'Blowdown';
pressurant.B = 4;
pressurant.pressure = 30e5 ; 

thruster.mass = 1.01 ;
thruster.number = 4 ;
thruster.chamber_pressure = 3.25e5 ;
thruster.chamber_pressure_min = 1.8e5 ; % [ Pa ] - Minimum combustion chamber pressure

tank.ratio = 1 ;
tank.number = 1 ;

[ sizing_SOSL_LANDING, M_PS_SOSL_LANDING ] = monopropellant_sizing( propellant, pressurant, tank, thruster )

% MONOPROPELLANT - SK - MR106E 22 N
% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.

% Orbiter + lander
dV_sk = 185*2;                             % Additional MAR-DV-020
Isp_sk = 235; 
m_dry_sk_orb_lan = 1279 ;                            % Sampling orbiter
m_prop_sk_orb_lan = preliminary_prop_mass(dV_sk,m_dry_sk_orb_lan,Isp_sk);

% Only orbiter
dV_sk = 185*2;                             % Additional MAR-DV-020
Isp_sk = 235; 
m_dry_sk_orb = 533;                            % Sampling orbiter
m_prop_sk_orb = preliminary_prop_mass(dV_sk,m_dry_sk_orb,Isp_sk);

% 
m_prop_sk = m_prop_sk_orb_lan + m_prop_sk_orb ;
%m_prop_sk = 205.8;

propellant.mass = m_prop_sk ;
propellant.rho = 1e3 ;
propellant.type = 'Cylindrical' ;

pressurant.feed_system = 'Blowdown';
pressurant.B = 4;
pressurant.pressure = 50e5 ;

thruster.mass = 0.635 ;
thruster.number = 16 ;
thruster.chamber_pressure = 18e5 ;
thruster.chamber_pressure_min = 4.5e5 ; % [ Pa ] - Minimum combustion pressure chamber

tank.ratio = 1 ;
tank.number = 1 ;

[ sizing_SOSL_SK, M_PS_SOSL_SK ] = monopropellant_sizing( propellant, pressurant, tank, thruster )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NON SAMPLING ORBITER - SAMPLING LANDER
close all, clear all, clc

%%%%%%%%
% GENERAL DATA
pressurant.R = 2077.3 ;
pressurant.gamma = 1.667 ;

tank.rho = 2780 ;
tank.sigma = 950e6 ;
%%%%%%%%

% BIPROPELLANT - 400N BI-PROPELLANT APOGEE THRUSTER
dV_int = 30 + 1.05*0.5e+3 + 1.05*0.9e+3 + 2*1.1e+3 + 10;  % plus MAR-DV-010 / MAR-DV-020 / MAR-DV-080 / MAR-DV-100
Isp_int = 321; m_dry_int = 928;                           % Orbiter + lander
m_prop_int = preliminary_prop_mass(dV_int,m_dry_int,Isp_int);
%m_prop_int = 2691.9;                    % from preliminary sizing

oxidizer.rho = 1443 ; 
oxidizer.OFratio = 1.65 ;
oxidizer.prop_mass = m_prop_int ;
oxidizer.type = 'Cylindrical' ;

fuel.rho = 900 ;

pressurant.feed_system = 'PressureRegulated' ;
pressurant.pressure = 300e5 ; 
pressurant.type = 'Spherical' ;

tank.ratio = 1 ;
tank.number = 3 ;

thruster.mass = 4.30 ;
thruster.number = 2 ;
thruster.chamber_pressure = 10.35e5 ;

[ M_PS_NSOSL_INT, sizing_ox_NSOSL_INT, sizing_fu_NSOSL_INT, sizing_gas_NSOSL_INT ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )

% MONOPROPELLANT - LANDING - MR107T 110 N
dV_land = 258*2;                            % Additional MAR-DV-020
Isp_land = 228; m_dry_land = 746;           % Sampling lander
m_prop_land = preliminary_prop_mass(dV_land,m_dry_land,Isp_land);
%m_prop_land = 95.33;

propellant.mass = m_prop_land ;
propellant.rho = 1e3 ;
propellant.type = 'Cylindrical' ;

pressurant.feed_system = 'Blowdown';
pressurant.B = 4;
pressurant.pressure = 30e5 ; 

thruster.mass = 1.01 ;
thruster.number = 4 ;
thruster.chamber_pressure = 3.25e5 ;
thruster.chamber_pressure_min = 1.8e5 ; % [ Pa ] - Minimum combustion chamber pressure

tank.ratio = 1 ;
tank.number = 1 ;

[ sizing_NSOSL_LANDING, M_PS_NSOSL_LANDING ] = monopropellant_sizing( propellant, pressurant, tank, thruster )

% MONOPROPELLANT - SK - MR111C 4 N
% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.

% Orbiter + lander
dV_sk = 275*2;                             % Additional MAR-DV-020
Isp_sk = 229; 
m_dry_sk_orb_lan = 928 ;                            % Sampling orbiter
m_prop_sk_orb_lan = preliminary_prop_mass(dV_sk,m_dry_sk_orb_lan,Isp_sk);

% Only orbiter
dV_sk = 275*2;                             % Additional MAR-DV-020
Isp_sk = 229; 
m_dry_sk_orb = 182;                            % Sampling orbiter
m_prop_sk_orb = preliminary_prop_mass(dV_sk,m_dry_sk_orb,Isp_sk);

% 
m_prop_sk = m_prop_sk_orb_lan + m_prop_sk_orb ;



% dV_sk = 550*2;                             % Additional MAR-DV-020
% Isp_sk = 229; m_dry_sk = 928;              % NS Orbiter + S Lander (for know, then change)
% m_prop_sk = preliminary_prop_mass(dV_sk,m_dry_sk,Isp_sk);
%m_prop_sk = 117.3 ;

propellant.mass = m_prop_sk ;
propellant.rho = 1e3 ;
propellant.type = 'Cylindrical' ;

pressurant.feed_system = 'Blowdown';
pressurant.B = 4;
pressurant.pressure = 50e5 ;

thruster.mass = 0.33 ;
thruster.number = 16 ;
thruster.chamber_pressure = 16e5 ;
thruster.chamber_pressure_min = 3.4e5 ; % [ Pa ] - Minimum combustion pressure chamber

tank.ratio = 1 ;
tank.number = 1 ;

[ sizing_NSOSL_SK, M_PS_NSOSL_SK ] = monopropellant_sizing( propellant, pressurant, tank, thruster )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%