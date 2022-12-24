clear all, close all, clc

% Compute the propellant mass for both modules of the selected arhitecture
% (SO-SL). Tanks are sized accordingly.

%%%%% Margins
n_flyby = 3 ;

MAR_010 = 1.05;          % Deterministic maneuvers
MAR_020 = 2;             % Stochastic
MAR_050 = 1.2;           % Apply on LANDING
MAR_030 = 2;             % Attitude maneuvres
MAR_080 = 30;            % [m/s] Launcher dispersion - apply only on primary
MAR_090 = 15*n_flyby;    % [m/s] For each flyby
MAR_100 = 10;            % [m/s] For Moon approach navigation maneuvre
MAR_hazard = 80;         % [kg] Sum to LANDING propellant mass
%%%%%

%%%%% % GENERAL DATA
% Helium pressurant
pressurant.R = 2077.3 ;
pressurant.gamma = 1.667 ;

% Titanium tanks
tank.rho = 2780 ;
tank.sigma = 950e6 ;
%%%%%



%% DATA

% Delta v
dV_prim = MAR_080 + MAR_010 * 1.0332e+3 + MAR_020 * 1.1e+3 + MAR_100 + MAR_090 ;
dV_sec = 370 * MAR_020 + 30 * MAR_030 ; 
dV_sec_c = dV_sec / 9 ; % SK and attitude when only orbiter
dV_sec_orb = dV_sec * 8 / 9 ; % SK and attitude when clamped configuration
dV_disp = 250 * MAR_020 ;

dV_lan = 244.6 * MAR_050 ;
dV_att_lan = 10 * MAR_030 ;

% Specific impulses
Isp_prim = 321 ; % Main thruster orbiter (400N bipropellant apogee motor)
Isp_sec = 292 ; % Secondary thruster orbiter (10N bipropellant)

Isp_lan = 228 ; % Main thruster lander (MD107T monopropellant)
Isp_att = 224 ; % Secondary thruster lander (1N monopropellant)

% Dry masses
m_dry_orb = 1468.5 ;
m_dry_lan = 495.3 + 20 ; % 20 kg for attitude hardware

%% OVERALL PROPELLANT MASS - both modules
[ m_prop_prim, m_prop_sec, m_prop_lan ] = propellant_mass( Isp_prim, dV_prim, dV_disp, ...
    Isp_sec, dV_sec_c, dV_sec_orb, m_dry_orb, Isp_lan, Isp_att, dV_lan, dV_att_lan, m_dry_lan, MAR_hazard ) ;

fprintf( 'Orbiter propellant mass: %f kg \n', ( m_prop_prim + m_prop_sec ) )
fprintf( 'Lander propellant mass: %f kg \n', ( m_prop_lan ) )

%% TANK SIZING - ORBITER

%% TANK SIZING - LANDER

