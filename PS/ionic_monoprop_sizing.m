% Architectures sizing considering ionic monopropellants for lander

%% SAMPLING ORBITER - SAMPLING LANDER
close all, clear all, clc

%%%%%%%%
% GENERAL DATA
pressurant.R = 2077.3 ;
pressurant.gamma = 1.667 ;

tank.rho = 2780 ;
tank.sigma = 950e6 ;

n_flyby = 4 ;

% MARGINS
MAR_010 = 1.05 ;                                                                  % Deterministic maneuver (only capture at Saturn)
MAR_020 = 2 ;                                                                     % Stocastic maneuver (only if not computed with statistical approach considering 3sigma )
MAR_030 = 2 ;                                                                     % Attitude control maneuvers (DO NOT CONSIDER FOR NOW)
MAR_050 = 1.2 ;                                                                   % Landing maneuver
MAR_080 = 30 ;                                                                    % [ m / s ] - Launcher dispersion (only interpalnetary)
MAR_090 = 15 * n_flyby ;                                                          % [ m / s ] - Gravity assist maneuver (15% for each flyby)
MAR_100 = 10 ;                                                                    % [ m / s ] - Moon approach navigation
MAR_hazard = 80 ;                                                                 % [ kg ] - Extra propellant to account for hazard maneuvers during landing 
%%%%%%%%


% BIPROPELLANT - 400N BI-PROPELLANT APOGEE THRUSTER

% Engine properties
Tmax = 440 ;                                                                      % [ N ] - Maximum thrust delivered (datasheet)
Isp_int = 321;                                                                    % [ s ] - Specific impulse (datasheet)
mdot_bip = 0.135 ;                                                                % [ kg / s ] Mass flow rate (datasheet)
m_dry_int = 1279;                                                                 % [ kg ] Orbiter + lander

%%%%%%%% Earth-Saturn transfer
dV_int_leg_nomar = 0.5e+3 ;

%n_maneuver = 4 ;
% dV_mis_int_leg =dV_int_leg_nomar * MAR_010 / n_maneuver ;

dV_mis_int_leg = ( MAR_010 * 0.9e3 - dV_int_leg_nomar ) * 0.01 ;                   % [ m / s ] - Maximum misalignment loss acceptable for interplanetary leg

dV_mis_int_leg = 1 ; % [ m / s ] - Lo imponiamo noi

m_prop_int_leg = preliminary_prop_mass( dV_int_leg_nomar * MAR_020, m_dry_int, Isp_int ) ; % [ kg ] - Propellant mass
m_int_leg = m_dry_int + m_prop_int_leg ;                                          % [ kg ] - Total mass during interplanetary leg

dt_int_leg = m_prop_int_leg / mdot_bip ;                                          % [ s ] - Maneuvering time

alpha_int_leg = acos( 1 - ( dV_mis_int_leg * m_int_leg ) / ( dt_int_leg * Tmax ) ) ; % Misalignment angle during interplanetary leg
alpha_int_leg = rad2deg( alpha_int_leg )

%%%%%%%% Capture at Saturn
dV_Sat_nomar = 0.9e3 ;
dV_mis_Sat = ( MAR_010 * 0.9e3 - dV_Sat_nomar ) * 0.01 ;                           % [ m / s ] - Maximum misalignment loss acceptable for capture at saturn

dV_mis_Sat = 1 ;

m_prop_Sat = preliminary_prop_mass( dV_Sat_nomar * MAR_010, m_dry_int, Isp_int);  % [ kg ] - Propellant mass
m_Sat = m_dry_int + m_prop_Sat ;                                                  % [ kg ] - Total mass during capture at Saturn

dt_Sat = m_prop_Sat / mdot_bip ;                                                  % [ s ] - Maneuvering time

alpha_Sat = acos( 1 - ( dV_mis_Sat * m_Sat ) / ( dt_Sat * Tmax ) ) ;              % [ rad ] - Misalignment angle during capture at Saturn
alpha_Sat = rad2deg( alpha_Sat )

%%%%%%%% Moon tour and orbit insertion
dV_EOI_nomar = 1.1e3 ;                                                            % [ m / s ] - Total impulse for moon tour and Encelaudt orbit insertion with no margins
dV_mis_EOI = ( 1.1e3 * MAR_020 - dV_EOI_nomar ) * 0.01 ;                           % [ m / s ] - Maximum misalignment loss acceptable for moon tour and EOI

dV_mis_EOI = 1 ;

m_prop_EOI = preliminary_prop_mass( dV_EOI_nomar * MAR_020, m_dry_int, Isp_int);  % [ kg ]
m_EOI = m_dry_int + m_prop_EOI ;                                                  % [ kg ] - Total mass during moon tour and Enceladus orbit insertion

dt_EOI = m_prop_EOI / mdot_bip ;                                                  % [ s ] -  Maneuvering time

alpha_EOI = acos( 1 - ( dV_mis_EOI * m_EOI ) / ( dt_EOI * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_EOI = rad2deg( alpha_EOI )

%%%%%%%% Cumulative bipropellant sizing
dV_int = MAR_020 * 0.5e+3 + MAR_010 * 0.9e3 + MAR_020 * 1.1e3 + MAR_080 + MAR_090 + MAR_100 ;  % [ m / s ] - Total cumulative impulse plus MAR-DV-010 / MAR-DV-020 / MAR-DV-080 / MAR-DV-100
m_prop_int = preliminary_prop_mass( dV_int, m_dry_int, Isp_int);                  % [ kg ] - Cumulative mass

oxidizer.rho = 1443 ;                                                             % [ kg / m ^ 3 ] - Oxidizer density
oxidizer.OFratio = 1.65 ;                                                         % Oxidizer-to-fuel ratio
oxidizer.prop_mass = m_prop_int ;
oxidizer.type = 'Cylindrical' ;

fuel.rho = 900 ;                                                                  % [ kg / m ^ 3 ] - Fuel density

pressurant.feed_system = 'PressureRegulated' ;
pressurant.pressure = 300e5 ; 
pressurant.type = 'Spherical' ;

tank.ratio = 1 ;
tank.number = 3 ;

thruster.mass = 4.30 ;
thruster.number = 2 ;
thruster.chamber_pressure = 10.35e5 ;

[ M_PS_SOSL_INT, sizing_ox_SOSL_INT, sizing_fu_SOSL_INT, sizing_gas_SOSL_INT ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster ) ;

% MONOPROPELLANT - LANDING - Conceptual idea (ADN-based thruster)
dV_land_nomar = 258 ;                          
dV_mis_land = ( 258 * MAR_050 - dV_land_nomar ) * 0.01 ;

dV_mis_land = 1 ;

Isp_land = 253; % [ s ] - LMP-103S propellant (ADN)
Tmax = 125 ; % [ N ] - Maximum thrust
mdot_land = 0.03925 ;

m_dry_land = 746;           % Sampling lander
m_prop_land = preliminary_prop_mass( dV_land_nomar * MAR_050, m_dry_land, Isp_land );
m_prop_land = m_prop_land + MAR_hazard ; % Accounts also for hazard maneuevres
m_land = m_dry_land + m_prop_land ;

dt_land = m_prop_land / mdot_land ;

alpha_land = acos( 1 - ( dV_mis_land * m_land ) / ( dt_land * Tmax ) ) ; % Misalignment angle during capture at Saturn
alpha_land = rad2deg( alpha_land )

propellant.mass = m_prop_land ;
propellant.rho = 1.24e3 ;
propellant.type = 'Cylindrical' ;

pressurant.feed_system = 'Blowdown';
pressurant.B = 4;
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

% MONOPROPELLANT - SK - MR106E 22 N
% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.

% Orbiter + lander
dV_sk_nomar = 185 ;                             % Additional MAR-DV-020
dV_mis_sk = ( 185 * MAR_020 - dV_sk_nomar ) * 0.01 ;

dV_mis_sk = 1 ;

Isp_sk = 235; 
Tmax = 30.7 ; % [ N ] - Maximum thrust (datasheet)
mdot_sk = 0.00905 ; % [kg/s]

m_dry_sk_orb_land = 1279 ;                            % Sampling orbiter
m_prop_sk_orb_land = preliminary_prop_mass( dV_sk_nomar * MAR_020 , m_dry_sk_orb_land, Isp_sk ) ;
m_sk_orb_land = m_dry_sk_orb_land + m_prop_sk_orb_land ;

dt_sk_orb_land = m_prop_sk_orb_land / mdot_sk ;

alpha_sk_orb_land = acos( 1 - ( dV_mis_sk * m_sk_orb_land ) / ( dt_sk_orb_land * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_land = rad2deg( alpha_sk_orb_land )

% Only orbiter
m_dry_sk_orb = 533;                            % Sampling orbiter
m_prop_sk_orb = preliminary_prop_mass( dV_sk_nomar * MAR_020, m_dry_sk_orb, Isp_sk ) ;
m_sk_orb = m_dry_sk_orb + m_prop_sk_orb ;

dt_sk_orb = m_prop_sk_orb / mdot_sk ;

alpha_sk_orb = acos( 1 - ( dV_mis_sk * m_sk_orb ) / ( dt_sk_orb * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb = rad2deg( alpha_sk_orb )

%
m_prop_sk = m_prop_sk_orb_land + m_prop_sk_orb ;
dt_sk_tot = dt_sk_orb + dt_sk_orb_land ;

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

n_flyby = 4 ;
%%%%%%%%

% MARGINS
MAR_010 = 1.05 ; % Deterministic maneuver (only capture at Saturn)
MAR_020 = 2 ; % Stocastic maneuver (only if not computed with statistical approach considering 3sigma )
MAR_030 = 2 ; % Attitude control maneuvers (DO NOT CONSIDER FOR NOW)
MAR_050 = 1.2 ; % Landing maneuver
MAR_080 = 30 ; % [ m / s ] - Launcher dispersion (only interpalnetary)
MAR_090 = 15 * n_flyby ; % [ m / s ] - Gravity assist maneuver (15% for each flyby)
MAR_100 = 10 ; % [ m / s ] - Moon approach navigation
MAR_hazard = 80 ; % [ kg ] - Extra propellant to account for hazard maneuvers during landing 


% BIPROPELLANT - 400N BI-PROPELLANT APOGEE THRUSTER

% Engine properties
Tmax = 440 ;                                                                      % [ N ] - Maximum thrust delivered (datasheet)
Isp_int = 321;                                                                    % [ s ] - Specific impulse (datasheet)
mdot_bip = 0.135 ;                                                                % [ kg / s ] Mass flow rate (datasheet)
m_dry_int = 928 ;                                                                 % [ kg ] Orbiter + lander

%%%%%%%% Earth-Saturn transfer
dV_int_leg_nomar = 0.5e+3 ;

%n_maneuver = 4 ;
% dV_mis_int_leg =dV_int_leg_nomar * MAR_010 / n_maneuver ;

dV_mis_int_leg = 1 ; % [ m / s ] - Lo imponiamo noi

m_prop_int_leg = preliminary_prop_mass( dV_int_leg_nomar * MAR_020, m_dry_int, Isp_int ) ; % [ kg ] - Propellant mass
m_int_leg = m_dry_int + m_prop_int_leg ;                                          % [ kg ] - Total mass during interplanetary leg

dt_int_leg = m_prop_int_leg / mdot_bip ;                                          % [ s ] - Maneuvering time

alpha_int_leg = acos( 1 - ( dV_mis_int_leg * m_int_leg ) / ( dt_int_leg * Tmax ) ) ; % Misalignment angle during interplanetary leg
alpha_int_leg = rad2deg( alpha_int_leg )

%%%%%%%% Capture at Saturn
dV_Sat_nomar = 0.9e3 ;
dV_mis_Sat = 1 ;

m_prop_Sat = preliminary_prop_mass( dV_Sat_nomar * MAR_010, m_dry_int, Isp_int);  % [ kg ] - Propellant mass
m_Sat = m_dry_int + m_prop_Sat ;                                                  % [ kg ] - Total mass during capture at Saturn

dt_Sat = m_prop_Sat / mdot_bip ;                                                  % [ s ] - Maneuvering time

alpha_Sat = acos( 1 - ( dV_mis_Sat * m_Sat ) / ( dt_Sat * Tmax ) ) ;              % [ rad ] - Misalignment angle during capture at Saturn
alpha_Sat = rad2deg( alpha_Sat )

%%%%%%%% Moon tour and orbit insertion
dV_EOI_nomar = 1.1e3 ;                                                            % [ m / s ] - Total impulse for moon tour and Encelaudt orbit insertion with no margins
dV_mis_EOI = 1 ;

m_prop_EOI = preliminary_prop_mass( dV_EOI_nomar * MAR_020, m_dry_int, Isp_int);  % [ kg ]
m_EOI = m_dry_int + m_prop_EOI ;                                                  % [ kg ] - Total mass during moon tour and Enceladus orbit insertion

dt_EOI = m_prop_EOI / mdot_bip ;                                                  % [ s ] -  Maneuvering time

alpha_EOI = acos( 1 - ( dV_mis_EOI * m_EOI ) / ( dt_EOI * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_EOI = rad2deg( alpha_EOI )

%%%%%%%% Cumulative bipropellant sizing
dV_int = MAR_020 * 0.5e+3 + MAR_010 * 0.9e3 + MAR_020 * 1.1e3 + MAR_080 + MAR_090 + MAR_100 ;  % [ m / s ] - Total cumulative impulse plus MAR-DV-010 / MAR-DV-020 / MAR-DV-080 / MAR-DV-100
m_prop_int = preliminary_prop_mass( dV_int, m_dry_int, Isp_int);                  % [ kg ] - Cumulative mass

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

[ M_PS_NSOSL_INT, sizing_ox_NSOSL_INT, sizing_fu_NSOSL_INT, sizing_gas_NSOSL_INT ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster ) ;

% MONOPROPELLANT - LANDING - MR107T 110 N
dV_land_nomar = 258 ;                          

dV_mis_land = 1 ;

Isp_land = 253; % [ s ] - LMP-103S propellant (ADN)
Tmax = 125 ; % [ N ] - Maximum thrust
mdot_land = 0.03925 ;

m_dry_land = 746;           % Sampling lander
m_prop_land = preliminary_prop_mass( dV_land_nomar * MAR_050, m_dry_land, Isp_land );
m_prop_land = m_prop_land + MAR_hazard ; % Accounts also for hazard maneuevres
m_land = m_dry_land + m_prop_land ;

dt_land = m_prop_land / mdot_land ;

alpha_land = acos( 1 - ( dV_mis_land * m_land ) / ( dt_land * Tmax ) ) ; % Misalignment angle during capture at Saturn
alpha_land = rad2deg( alpha_land )

%
propellant.mass = m_prop_land ;
propellant.rho = 1.24e3 ;
propellant.type = 'Cylindrical' ;

pressurant.feed_system = 'Blowdown';
pressurant.B = 4;
pressurant.pressure = 50e5 ; 

thruster.mass = 1.01 ;
thruster.number = 4 ;
thruster.chamber_pressure = 3.25e5 ;
thruster.chamber_pressure_min = 1.8e5 ; % [ Pa ] - Minimum combustion chamber pressure

tank.ratio = 1 ;
tank.number = 1 ;

[ sizing_NSOSL_LANDING, M_PS_NSOSL_LANDING ] = monopropellant_sizing( propellant, pressurant, tank, thruster ) ;

% MONOPROPELLANT - SK - MR111C 4 N
% Assumption: half of the deltav for SK is used when orbiter and lander are
% clamped, while the other half is used only by the orbiter.

% Orbiter + lander
dV_sk_nomar = 275 ;                             % Additional MAR-DV-020
dV_mis_sk = 1 ;

Isp_sk = 229; 
mdot_sk = 0.9 ;
m_dry_sk_orb_land = 928 ;                            % Sampling orbiter
m_prop_sk_orb_land = preliminary_prop_mass( dV_sk_nomar * MAR_020 , m_dry_sk_orb_land, Isp_sk ) ;
m_sk_orb_land = m_dry_sk_orb_land + m_prop_sk_orb_land ;

dt_sk_orb_land = m_prop_sk_orb_land / mdot_sk ;

alpha_sk_orb_land = acos( 1 - ( dV_mis_sk * m_sk_orb_land ) / ( dt_sk_orb_land * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb_land = rad2deg( alpha_sk_orb_land )

% Only orbiter
m_dry_sk_orb = 182 ;                            % Sampling orbiter

m_prop_sk_orb = preliminary_prop_mass( dV_sk_nomar * MAR_020, m_dry_sk_orb, Isp_sk ) ;
m_sk_orb = m_dry_sk_orb + m_prop_sk_orb ;

dt_sk_orb = m_prop_sk_orb / mdot_sk ;

alpha_sk_orb = acos( 1 - ( dV_mis_sk * m_sk_orb ) / ( dt_sk_orb * Tmax ) ) ;              % [ rad ] - Misalignment angle during interplanetary leg
alpha_sk_orb = rad2deg( alpha_sk_orb )

% 
m_prop_sk = m_prop_sk_orb_land + m_prop_sk_orb ;


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