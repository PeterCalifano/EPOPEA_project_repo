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

n_flyby = 4;             % Number of flybys

% Margins
MAR_010 = 1.05;          % Deterministic maneuvers
MAR_020 = 2;             % Stochastic
MAR_050 = 1.2;           % Apply on LANDING
MAR_030 = 2;             % Attitude maneuvres
MAR_080 = 30;            % [m/s] Launcher dispersion - apply only on primary
MAR_090 = 15*n_flyby;    % [m/s] For each flyby
MAR_100 = 10;            % [m/s] For Moon approach navigation maneuvre
MAR_hazard = 80;         % [kg] Sum to LANDING propellant mass

% ----------- BIPROPELLANT - 400N BI-PROPELLANT APOGEE THRUSTER -----------

dV_int = MAR_080 + MAR_010*2.25e+3 + MAR_020*1.1e+3 + MAR_100 + MAR_090;   % [m/s] - Interplanetary + capture + moon toru and EOI
Isp_int = 321; 
Tmax = 440 ; % [N] - Maximum thrust
m_dry_int = 1279;                                                          % Orbiter + lander
mdot_bip = 0.135 ;                                                         % [ kg / s ] Mass flow rate (datasheet)
m_prop_int = preliminary_prop_mass(dV_int,m_dry_int,Isp_int);

% Compute maneuvering  time
m_prop_int1 = preliminary_prop_mass(MAR_010*1.79e3,m_dry_int,Isp_int) ;
m_int = m_prop_int1 + m_dry_int ;
dt_int = m_prop_int1 / mdot_bip ;

m_prop_cap = preliminary_prop_mass(MAR_010*0.46e3,m_dry_int,Isp_int) ;
m_cap = m_prop_cap + m_dry_int ;
dt_cap = m_prop_cap / mdot_bip ;

m_prop_EOI = preliminary_prop_mass(MAR_020*1.1e+3,m_dry_int,Isp_int) ;
m_EOI = m_prop_EOI + m_dry_int ;
dt_EOI = m_prop_EOI / mdot_bip ;

% Compute misalignment angle
dV_mis = 1 ; % [m/s] - Assumption

alpha_int_leg_SOSL = acos( 1 - ( dV_mis * m_int ) / ( dt_int * Tmax ) ) ; % Misalignment angle during interplanetary leg
alpha_int_leg_SOSL = rad2deg( alpha_int_leg_SOSL ) ;

alpha_cap_SOSL = acos( 1 - ( dV_mis * m_cap ) / ( dt_cap * Tmax ) ) ;              % [ rad ] - Misalignment angle during capture at Saturn
alpha_cap_SOSL = rad2deg( alpha_cap_SOSL ) ;

alpha_EOI_SOSL = acos( 1 - ( dV_mis * m_EOI ) / ( dt_EOI * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_EOI_SOSL = rad2deg( alpha_EOI_SOSL ) ;

% Bipropellant sizing
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

[ M_PS_SOSL_INT, sizing_ox_SOSL_INT, sizing_fu_SOSL_INT, sizing_gas_SOSL_INT ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster ) ;

% ---------------- MONOPROPELLANT - LANDING - MR107T 110 N ----------------
dV_land = 258*MAR_050;                                                     % Additional MAR-DV-020
Isp_land = 228; 
Tmax = 125 ;
m_dry_land = 746;                                                          % Sampling lander
mdot_land = 0.03925;                                                       % [kg/s]
% m_prop_land = preliminary_prop_mass(dV_land,m_dry_land,Isp_land) + MAR_hazard;
m_prop_land = 65 * MAR_020 + MAR_hazard ; % From MA 
m_land = m_prop_land + m_dry_land ;

% Compute maneuvering time
dt_land = m_prop_land / mdot_land;

% Compute misalignment angles
dV_mis = 1 ; % [m/s] - Assumption
alpha_land = acos( 1 - ( dV_mis * m_land ) / ( dt_land * Tmax ) ) ; % Misalignment angle during capture at Saturn
alpha_land = rad2deg( alpha_land ) ;

% Monopropellant sizing
propellant.mass = m_prop_land ;
propellant.rho = 1e3 ;
propellant.type = 'Cylindrical' ;

pressurant.feed_system = 'Blowdown';
pressurant.B = 4;
pressurant.pressure = 30e5 ; 

thruster.mass = 1.01 ;
thruster.number = 4 ;
thruster.chamber_pressure = 3.25e5 ;
thruster.chamber_pressure_min = 1.8e5 ;                                    % [ Pa ] - Minimum combustion chamber pressure

tank.ratio = 1 ;
tank.number = 1 ;

[ sizing_SOSL_LANDING, M_PS_SOSL_LANDING ] = monopropellant_sizing( propellant, pressurant, tank, thruster ) ;

% ------------------- MONOPROPELLANT - SK - MR106E 22 N -------------------

% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.

% Orbiter + lander
dV_sk = 185*MAR_020;                             % Additional MAR-DV-020
Isp_sk = 235; 
Tmax = 30.7 ;
m_dry_sk_orb_lan = 1279 ;                            % Sampling orbiter
mdot_sk = 0.00905 ;
m_prop_sk_orb_lan = preliminary_prop_mass(dV_sk,m_dry_sk_orb_lan,Isp_sk);
m_sk_orb_lan = m_prop_sk_orb_lan + m_dry_sk_orb_lan ;

% Only orbiter
dV_sk = 185*MAR_020;                             % Additional MAR-DV-020
Isp_sk = 235; 
m_dry_sk_orb = 533;                            % Sampling orbiter
m_prop_sk_orb = preliminary_prop_mass(dV_sk,m_dry_sk_orb,Isp_sk);
m_sk_orb = m_prop_sk_orb + m_dry_sk_orb ;

% Total
m_prop_sk = m_prop_sk_orb_lan + m_prop_sk_orb ;

% Compute maneuvering time
dt_sk1 = m_prop_sk_orb_lan / mdot_sk;
dt_sk2 = m_prop_sk_orb / mdot_sk ; 

% Compute misalignment angle
dV_mis = 1 ; % [m/s] - Assumption
alpha_sk_orb_land_SOSL = acos( 1 - ( dV_mis * m_sk_orb_lan ) / ( dt_sk1 * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_land_SOSL = rad2deg( alpha_sk_orb_land_SOSL ) ;

alpha_sk_orb_SOSL = acos( 1 - ( dV_mis * m_sk_orb ) / ( dt_sk2 * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_SOSL = rad2deg( alpha_sk_orb_SOSL ) ;

% Monopropellant sizing
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

[ sizing_SOSL_SK, M_PS_SOSL_SK ] = monopropellant_sizing( propellant, pressurant, tank, thruster ) ;

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

n_flyby = 4;             % Number of flybys

% Margins
MAR_010 = 1.05;          % Deterministic maneuvers
MAR_020 = 2;             % Stochastic
MAR_050 = 1.2;           % Apply on LANDING
MAR_030 = 2;             % Attitude maneuvres
MAR_080 = 30;            % [m/s] Launcher dispersion - apply only on primary
MAR_090 = 15*n_flyby;    % [m/s] For each flyby
MAR_100 = 10;            % [m/s] For Moon approach navigation maneuvre
MAR_hazard = 80;         % [kg] Sum to LANDING propellant mass

% ----------- BIPROPELLANT - 400N BI-PROPELLANT APOGEE THRUSTER ----------- 

dV_int = MAR_080 + MAR_010*2.25e+3 + MAR_020*1.1e+3 + MAR_100 + MAR_090;  % plus MAR-DV-010 / MAR-DV-020 / MAR-DV-080 / MAR-DV-100
Isp_int = 321;
Tmax = 440;
m_dry_int = 928;                                                           % Orbiter + lander
mdot_bip = 0.135 ;
m_prop_int = preliminary_prop_mass(dV_int,m_dry_int,Isp_int);

% Compute maneuvering time
m_prop_int1 = preliminary_prop_mass(MAR_010*1.79e+3,m_dry_int,Isp_int) ;
m_int = m_prop_int1 + m_dry_int ;
dt_int = m_prop_int1 / mdot_bip ;

m_prop_cap = preliminary_prop_mass(MAR_010*0.46e+3,m_dry_int,Isp_int) ;
m_cap = m_prop_cap + m_dry_int ;
dt_cap = m_prop_cap / mdot_bip ;

m_prop_EOI = preliminary_prop_mass(MAR_020*1.1e+3,m_dry_int,Isp_int) ;
m_EOI = m_prop_EOI + m_dry_int ;
dt_EOI = m_prop_EOI / mdot_bip ;

% Compute misalignment angle
dV_mis = 1 ; % [m/s] - Assumption

alpha_int_leg_NSOSL = acos( 1 - ( dV_mis * m_int ) / ( dt_int * Tmax ) ) ; % Misalignment angle during interplanetary leg
alpha_int_leg_NSOSL = rad2deg( alpha_int_leg_NSOSL ) ;

alpha_cap_NSOSL = acos( 1 - ( dV_mis * m_cap ) / ( dt_cap * Tmax ) ) ;              % [ rad ] - Misalignment angle during capture at Saturn
alpha_cap_NSOSL = rad2deg( alpha_cap_NSOSL ) ;

alpha_EOI_NSOSL = acos( 1 - ( dV_mis * m_EOI ) / ( dt_EOI * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_EOI_NSOSL = rad2deg( alpha_EOI_NSOSL ) ;

% Bipropellant sizing
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

[ M_PS_NSOSL_INT, sizing_ox_NSOSL_INT, sizing_fu_NSOSL_INT, sizing_gas_NSOSL_INT ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster ) ;

% ---------------- MONOPROPELLANT - LANDING - MR107T 110 N ---------------- 
dV_land = 258*MAR_050;                            % Additional MAR-DV-020
Isp_land = 228; 
Tmax = 125 ;
m_dry_land = 746;           % Sampling lander
mdot_land = 0.03925;  % [kg/s]
% m_prop_land = preliminary_prop_mass(dV_land,m_dry_land,Isp_land) + MAR_hazard;
m_prop_land = 65 * MAR_020 + MAR_hazard ; % From MA 
m_land = m_prop_land + m_dry_land ;

% Compute maneuvering time
dt_land = m_prop_land / mdot_land;

% Compute misalignment angle
dV_mis = 1 ; % [m/s] - Assumption
alpha_land_NSOSL = acos( 1 - ( dV_mis * m_land ) / ( dt_land * Tmax ) ) ; % Misalignment angle during capture at Saturn
alpha_land_NSOSL = rad2deg( alpha_land_NSOSL ) ;

% Monopropellant sizing
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

[ sizing_NSOSL_LANDING, M_PS_NSOSL_LANDING ] = monopropellant_sizing( propellant, pressurant, tank, thruster ) ;

% ---------------- MONOPROPELLANT - SK - MR111C 4 N
% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.

% Orbiter + lander
dV_sk = 275*MAR_020;                             % Additional MAR-DV-020
Isp_sk = 229; 
m_dry_sk_orb_lan = 928 ;                            % Sampling orbiter
mdot_sk = 0.0009;
m_prop_sk_orb_lan = preliminary_prop_mass(dV_sk,m_dry_sk_orb_lan,Isp_sk);
m_sk_orb_lan = m_prop_sk_orb_lan + m_dry_sk_orb_lan ;

% Only orbiter
dV_sk = 275*MAR_020;                             % Additional MAR-DV-020
Isp_sk = 229; 
m_dry_sk_orb = 182;                            % Sampling orbiter
m_prop_sk_orb = preliminary_prop_mass(dV_sk,m_dry_sk_orb,Isp_sk);
m_sk_orb = m_prop_sk_orb + m_dry_sk_orb ;

m_prop_sk = m_prop_sk_orb_lan + m_prop_sk_orb ;

% Compute maneuvering time
dt_sk1 = m_prop_sk_orb_lan / mdot_sk;
dt_sk2 = m_prop_sk_orb / mdot_sk ; 

% Compute misalignment angle
dV_mis = 1 ; % [m/s] - Assumption
alpha_sk_orb_land_NSOSL = acos( 1 - ( dV_mis * m_sk_orb_lan ) / ( dt_sk1 * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_land_NSOSL = rad2deg( alpha_sk_orb_land_NSOSL ) ;

alpha_sk_orb_NSOSL = acos( 1 - ( dV_mis * m_sk_orb ) / ( dt_sk2 * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_NSOSL = rad2deg( alpha_sk_orb_NSOSL ) ;

% Monopropellant sizing
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

[ sizing_NSOSL_SK, M_PS_NSOSL_SK ] = monopropellant_sizing( propellant, pressurant, tank, thruster ) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%