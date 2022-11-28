% Architectures sizing considering a dual mode DM - pressure regulated 
% for the primary propulsion (MAIN ENGINE) and 
% blow-down for the secondary propulsion (SMALLER THRUSTERS)

% BIPROPELLANT thrusters:
% Main Engine = ARIANE 400N BI-PROPELLANT APOGEE THRUSTER
% Smaller Thrusters = ARIANE 10N BI-PROPELLANT THRUSTER

%% SAMPLING ORBITER - SAMPLING LANDER
close all, clear all, clc

%%%%%%%%
% GENERAL DATA
pressurant.R = 2077.3 ;
pressurant.gamma = 1.667 ;

tank.rho = 2780 ;
tank.sigma = 950e6 ;

mdot_reg = 0.135;
mdot_blow = 0.0035;
%%%%%%%%
n_flyby = 3;             % Number of flybys

% Margins
MAR_010 = 1.05;          % Deterministic maneuvers
MAR_020 = 2;             % Stochastic
MAR_050 = 1.2;           % Apply on LANDING
MAR_030 = 2;             % Attitude maneuvres
MAR_080 = 30;            % [m/s] Launcher dispersion - apply only on primary
MAR_090 = 15*n_flyby;    % [m/s] For each flyby
MAR_100 = 10;            % [m/s] For Moon approach navigation maneuvre
MAR_hazard = 80;         % [kg] Sum to LANDING propellant mass

% ----- PRESSURE REGULATED: ARIANE 400N BI-PROPELLANT APOGEE THRUSTER -----
dV_reg = MAR_080 + MAR_010*1.0158e+3 + MAR_020*1.1e+3 + MAR_100 + MAR_090;  % plus MAR-DV-010 / MAR-DV-020 / MAR-DV-080 / MAR-DV-100
Isp_reg = 321; 
Tmax = 440 ;
m_dry_int = 1279;                          % Orbiter + lander
m_prop_reg = preliminary_prop_mass(dV_reg,m_dry_int,Isp_reg);

% Compute maneuvering time
dV_int_leg_nomar = 1.0158 - 0.248 ;
m_prop_int1 = preliminary_prop_mass(MAR_010*dV_int_leg_nomar,m_dry_int,Isp_reg) ;
m_int = m_prop_int1 + m_dry_int ;
dt_int = m_prop_int1 / mdot_reg ;

dV_Sat_nomar = 0.248 ;
m_prop_cap = preliminary_prop_mass(MAR_010*dV_Sat_nomar,m_dry_int,Isp_reg) ;
m_cap = m_prop_cap + m_dry_int ;
dt_cap = m_prop_cap / mdot_reg ;

m_prop_EOI = preliminary_prop_mass(MAR_020*1.1e+3,m_dry_int,Isp_reg) ;
m_eoi = m_prop_EOI + m_dry_int ;
dt_EOI = m_prop_EOI / mdot_reg ;

% Compute misalignment angle
dV_mis = 1 ; % [m/s] - Assumption

alpha_int_leg_SOSL = acos( 1 - ( dV_mis * m_int ) / ( dt_int * Tmax ) ) ; % Misalignment angle during interplanetary leg
alpha_int_leg_SOSL = rad2deg( alpha_int_leg_SOSL ) ;

alpha_cap_SOSL = acos( 1 - ( dV_mis * m_cap ) / ( dt_cap * Tmax ) ) ;              % [ rad ] - Misalignment angle during capture at Saturn
alpha_cap_SOSL = rad2deg( alpha_cap_SOSL ) ;

alpha_EOI_SOSL = acos( 1 - ( dV_mis * m_eoi ) / ( dt_EOI * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_EOI_SOSL = rad2deg( alpha_EOI_SOSL ) ;

% -------------- BLOWDOWN: ARIANE 10N BI-PROPELLANT THRUSTER --------------

% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.

% Orbiter + lander
dV_blow = 185*MAR_020;                             % Additional MAR-DV-020
Isp_blow = 292; 
Tmax = 12.5 ;
m_dry_sk_orb_lan = 1279 ;                    % Sampling orbiter
m_prop_sk_orb_lan = preliminary_prop_mass(dV_blow,m_dry_sk_orb_lan,Isp_blow);
m_sk_orb_lan = m_prop_sk_orb_lan + m_dry_sk_orb_lan ;

% Only orbiter
dV_blow = 185*MAR_020;                             % Additional MAR-DV-020
Isp_blow = 292; 
m_dry_sk_orb = 533;                          % Sampling orbiter
m_prop_sk_orb = preliminary_prop_mass(dV_blow,m_dry_sk_orb,Isp_blow);
m_sk_orb = m_prop_sk_orb + m_dry_sk_orb ;

m_prop_blow = m_prop_sk_orb_lan + m_prop_sk_orb ;

% Compute maneuvering time
dt_sk1 = m_prop_sk_orb_lan / mdot_blow ;
dt_sk2 = m_prop_sk_orb / mdot_blow ;

% Compute misalignment angle
dV_mis = 1 ; % [m/s] - Assumption
alpha_sk_orb_land_SOSL = acos( 1 - ( dV_mis * m_sk_orb_lan ) / ( dt_sk1 * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_land_SOSL = rad2deg( alpha_sk_orb_land_SOSL ) ;

alpha_sk_orb_SOSL = acos( 1 - ( dV_mis * m_sk_orb ) / ( dt_sk2 * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_SOSL = rad2deg( alpha_sk_orb_SOSL ) ;

% DATA
oxidizer.rho = 1443 ; oxidizer.OFratio_reg = 1.65 ; oxidizer.OFratio_blow = 1.60 ;
oxidizer.prop_mass_reg = m_prop_reg ; oxidizer.prop_mass_blow = m_prop_blow ;
oxidizer.type = 'Cylindrical' ;
fuel.rho = 900 ;

pressurant.pressure = 300e5 ; pressurant.type = 'Spherical' ;
pressurant.B = 4;
tank.ratio = 1 ; tank.number = 3 ;

thruster.mass_ME = 4.30 ; thruster.number_ME = 2 ;
thruster.mass_ST = 0.65 ; thruster.number_ST = 16 ;
thruster.chamber_pressure_reg = 10.35e5 ; thruster.chamber_pressure_blow = 9e5 ;

disp('S ORBITER - S LANDER:')
[ M_PS_SOSL, sizing_ox_SOSL, sizing_fu_SOSL, sizing_gas_SOSL ] = bipropellant_DM( oxidizer, fuel, pressurant, tank, thruster ) ;

%% NON SAMPLING ORBITER - SAMPLING LANDER
close all, clear all, clc

%%%%%%%%
% GENERAL DATA
pressurant.R = 2077.3 ;
pressurant.gamma = 1.667 ;

tank.rho = 2780 ;
tank.sigma = 950e6 ;

mdot_reg = 0.135;
mdot_blow = 0.0035;
%%%%%%%%

n_flyby = 3 ;

% Margins
MAR_010 = 1.05;          % Deterministic maneuvers
MAR_020 = 2;             % Stochastic
MAR_050 = 1.2;           % Apply on LANDING
MAR_030 = 2;             % Attitude maneuvres
MAR_080 = 30;            % [m/s] Launcher dispersion - apply only on primary
MAR_090 = 15*n_flyby;    % [m/s] For each flyby
MAR_100 = 10;            % [m/s] For Moon approach navigation maneuvre
MAR_hazard = 80;         % [kg] Sum to LANDING propellant mass

% ----- PRESSURE REGULATED: ARIANE 400N BI-PROPELLANT APOGEE THRUSTER -----
dV_reg = MAR_080 + MAR_010 *1.0158e+3 + MAR_020*1.1e+3 + MAR_100 + MAR_090;  % plus MAR-DV-010 / MAR-DV-020 / MAR-DV-080 / MAR-DV-100
Isp_reg = 321; 
Tmax = 440;
m_dry_int = 928;                           % Orbiter + lander
m_prop_reg = preliminary_prop_mass(dV_reg,m_dry_int,Isp_reg);

% Compute maneuvering time
dV_int_leg_nomar = 1.0158 - 0.248 ;
m_prop_int = preliminary_prop_mass(MAR_010*dV_int_leg_nomar,m_dry_int,Isp_reg) ;
m_int = m_prop_int + m_dry_int ;
dt_int = m_prop_int / mdot_reg ;

dV_Sat_nomar = 0.248 ;
m_prop_cap = preliminary_prop_mass(MAR_010*dV_Sat_nomar,m_dry_int,Isp_reg) ;
m_cap = m_prop_cap + m_dry_int ;
dt_cap = m_prop_cap / mdot_reg ;

m_prop_eoi = preliminary_prop_mass(2*1.1e+3,m_dry_int,Isp_reg) ;
m_eoi = m_prop_eoi + m_dry_int ;
dt_eoi = m_prop_eoi / mdot_reg ;

% Compute misalignment angle
dV_mis = 1 ; % [m/s] - Assumption

alpha_int_leg_NSOSL = acos( 1 - ( dV_mis * m_int ) / ( dt_int * Tmax ) ) ; % Misalignment angle during interplanetary leg
alpha_int_leg_NSOSL = rad2deg( alpha_int_leg_NSOSL ) ;

alpha_cap_NSOSL = acos( 1 - ( dV_mis * m_cap ) / ( dt_cap * Tmax ) ) ;              % [ rad ] - Misalignment angle during capture at Saturn
alpha_cap_NSOSL = rad2deg( alpha_cap_NSOSL ) ;

alpha_EOI_NSOSL = acos( 1 - ( dV_mis * m_eoi ) / ( dt_eoi * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_EOI_NSOSL = rad2deg( alpha_EOI_NSOSL ) ;

% -------------- BLOWDOWN: ARIANE 10N BI-PROPELLANT THRUSTER --------------

% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.

% Orbiter + lander
dV_sk = 275*MAR_020;                             % Additional MAR-DV-020
Isp_sk = 292; 
Tmax = 12.5 ;
m_dry_sk_orb_lan = 928 ;                            % Sampling orbiter
m_prop_sk_orb_lan = preliminary_prop_mass(dV_sk,m_dry_sk_orb_lan,Isp_sk);
m_sk_orb_lan = m_prop_sk_orb_lan + m_dry_sk_orb_lan ;

% Only orbiter
dV_sk = 275*MAR_020;                             % Additional MAR-DV-020
Isp_sk = 292; 
m_dry_sk_orb = 182;                            % Sampling orbiter
m_prop_sk_orb = preliminary_prop_mass(dV_sk,m_dry_sk_orb,Isp_sk);
m_sk_orb = m_prop_sk_orb + m_dry_sk_orb ;

m_prop_blow = m_prop_sk_orb_lan + m_prop_sk_orb ;

% Compute maneuvering time
dt_sk1 = m_prop_sk_orb_lan / mdot_blow ;
dt_sk2 = m_prop_sk_orb / mdot_blow ; 

% Compute misalignment angle
dV_mis = 1 ; % [m/s] - Assumption
alpha_sk_orb_land_NSOSL = acos( 1 - ( dV_mis * m_sk_orb_lan ) / ( dt_sk1 * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_land_NSOSL = rad2deg( alpha_sk_orb_land_NSOSL ) ;

alpha_sk_orb_NSOSL = acos( 1 - ( dV_mis * m_sk_orb ) / ( dt_sk2 * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_NSOSL = rad2deg( alpha_sk_orb_NSOSL ) ;

% DATA
oxidizer.rho = 1443 ; oxidizer.OFratio_reg = 1.65 ; oxidizer.OFratio_blow = 1.60 ;
oxidizer.prop_mass_reg = m_prop_reg ; oxidizer.prop_mass_blow = m_prop_blow ;
oxidizer.type = 'Cylindrical' ;
fuel.rho = 900 ;

pressurant.pressure = 300e5 ; pressurant.type = 'Spherical' ;
pressurant.B = 4;
tank.ratio = 1 ; tank.number = 3 ;

thruster.mass_ME = 4.30 ; thruster.number_ME = 2 ;
thruster.mass_ST = 0.65 ; thruster.number_ST = 16 ;
thruster.chamber_pressure_reg = 10.35e5 ; thruster.chamber_pressure_blow = 9e5 ;

disp('NS ORBITER - S LANDER:')
[ M_PS_NSOSL, sizing_ox_NSOSL, sizing_fu_NSOSL, sizing_gas_NSOSL] = bipropellant_DM( oxidizer, fuel, pressurant, tank, thruster ) ;

