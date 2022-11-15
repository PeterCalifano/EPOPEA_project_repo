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
%%%%%%%%

% PRESSURE REGULATED: ARIANE 400N BI-PROPELLANT APOGEE THRUSTER
dV_reg = 30 + 1.05*0.5e+3 + 1.05*0.9e+3 + 2*1.1e+3 + 10;  % plus MAR-DV-010 / MAR-DV-020 / MAR-DV-080 / MAR-DV-100
Isp_reg = 321; m_dry_int = 1279;                          % Orbiter + lander
m_prop_reg = preliminary_prop_mass(dV_reg,m_dry_int,Isp_reg);

% BLOWDOWN: ARIANE 10N BI-PROPELLANT THRUSTER
% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.
% Orbiter + lander
dV_blow = 185*2;                             % Additional MAR-DV-020
Isp_blow = 292; 
m_dry_sk_orb_lan = 1279 ;                    % Sampling orbiter
m_prop_sk_orb_lan = preliminary_prop_mass(dV_blow,m_dry_sk_orb_lan,Isp_blow);
% Only orbiter
dV_blow = 185*2;                             % Additional MAR-DV-020
Isp_blow = 292; 
m_dry_sk_orb = 533;                          % Sampling orbiter
m_prop_sk_orb = preliminary_prop_mass(dV_blow,m_dry_sk_orb,Isp_blow);

m_prop_blow = m_prop_sk_orb_lan + m_prop_sk_orb ;

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
[ M_PS_SS, sizing_ox_SS, sizing_fu_SS, sizing_SS ] = bipropellant_DM( oxidizer, fuel, pressurant, tank, thruster )

%% NON SAMPLING ORBITER - SAMPLING LANDER
close all, clear all

%%%%%%%%
% GENERAL DATA
pressurant.R = 2077.3 ;
pressurant.gamma = 1.667 ;

tank.rho = 2780 ;
tank.sigma = 950e6 ;
%%%%%%%%

% PRESSURE REGULATED: ARIANE 400N BI-PROPELLANT APOGEE THRUSTER
dV_reg = 30 + 1.05*0.5e+3 + 1.05*0.9e+3 + 2*1.1e+3 + 10;  % plus MAR-DV-010 / MAR-DV-020 / MAR-DV-080 / MAR-DV-100
Isp_reg = 321; m_dry_int = 928;                           % Orbiter + lander
m_prop_reg = preliminary_prop_mass(dV_reg,m_dry_int,Isp_reg);

% BLOWDOWN: ARIANE 10N BI-PROPELLANT THRUSTER
% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.
% Orbiter + lander
dV_sk = 275*2;                             % Additional MAR-DV-020
Isp_sk = 292; 
m_dry_sk_orb_lan = 928 ;                            % Sampling orbiter
m_prop_sk_orb_lan = preliminary_prop_mass(dV_sk,m_dry_sk_orb_lan,Isp_sk);

% Only orbiter
dV_sk = 275*2;                             % Additional MAR-DV-020
Isp_sk = 292; 
m_dry_sk_orb = 182;                            % Sampling orbiter
m_prop_sk_orb = preliminary_prop_mass(dV_sk,m_dry_sk_orb,Isp_sk);

m_prop_blow = m_prop_sk_orb_lan + m_prop_sk_orb ;

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
[ M_PS_NS, sizing_ox_NS, sizing_fu_NS, sizing_NS ] = bipropellant_DM( oxidizer, fuel, pressurant, tank, thruster )

