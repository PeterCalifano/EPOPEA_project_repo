clear all, close all, clc

% LANDER SIZING
% IONIC MONOPROPELLANT (ADN based thruster)

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

%%%%%%%%
% GENERAL DATA
pressurant.R = 2077.3 ;
pressurant.gamma = 1.667 ;

tank.rho = 2780 ;
tank.sigma = 950e6 ;
%%%%%%%%

% Mass sizing
dV_land = 245.5 * MAR_050 ;

dV_mis_land = 1 ;

Isp_land = 253; % [ s ] - LMP-103S propellant (ADN)
Tmax = 125 ; % [ N ] - Maximum thrust
mdot_land = 0.03925 ;

m_dry_land = 452.16 ;          % Sampling lander
m_prop_land = preliminary_prop_mass( dV_land, m_dry_land, Isp_land );
m_prop_land = m_prop_land + MAR_hazard ; % Accounts also for hazard maneuevres
m_land = m_dry_land + m_prop_land ;

dt_land = m_prop_land / mdot_land ;

alpha_land = acos( 1 - ( dV_mis_land * m_land ) / ( dt_land * Tmax ) ) ; % Misalignment angle during capture at Saturn
alpha_land = rad2deg( alpha_land ) ;

propellant.mass = m_prop_land ;
propellant.rho = 1.24e3 ;
propellant.type = 'Cylindrical' ;

pressurant.feed_system = 'Blowdown' ;
pressurant.B = 4 ;
pressurant.pressure = 50e5 ;

% We don-'t have a reference thruster yet, therefore I leave the data of the 
% hydrazine monoprop thrusters
thruster.mass = 1.01 ;
thruster.number = 4 ;
thruster.chamber_pressure = 3.25e5 ;
thruster.chamber_pressure_min = 1.8e5 ; % [ Pa ] - Minimum combustion chamber pressure

tank.ratio = 1 ;
tank.number = 1 ;

[ sizing_SOSL_LANDING, M_PS_SOSL_LANDING ] = monopropellant_sizing( propellant, pressurant, tank, thruster ) ;
