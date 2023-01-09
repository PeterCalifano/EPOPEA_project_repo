clear all, close all, clc

% Compute the propellant mass for both modules of the selected arhitecture
% (SO-SL). Tanks are sized accordingly.

%%%%% Margins
n_flyby = 3 ;

MAR_010 = 1.05;          % Deterministic maneuvers
MAR_020 = 2;             % Stochastic
MAR_030 = 2;             % Attitude maneuvers
MAR_050 = 1.2;           % Apply on LANDING
MAR_080 = 30;            % [m/s] Launcher dispersion - apply only on primary
MAR_090 = 15*n_flyby;    % [m/s] For each flyby
MAR_100 = 10;            % [m/s] For Moon approach navigation maneuver
MAR_attitude = 1.5;      % Propellant mass for attitude
MAR_hazard = 80;         % [kg] Sum to LANDING propellant mass
%%%%%

                % %% DATA - STANDARD CONFIGURATION
                % 
                % % Lander
                % lander.dV_des = 244.6 * MAR_050 ; % [m/s] - Delta v for descending phase
                % lander.dV_att = 10 * MAR_030 ; % [m/s] - Delta v for attitude control during landing phase
                % lander.Isp_main = 228 ; % [s] - Specific impulse of main thrusters for landing
                % lander.Isp_att = 224 ; % [s] - Specific impulse of attitude thrusters for landing
                % 
                % % Orbiter
                % dV_sk = 370 * MAR_020 ; % [m/s] - Cumulative delta v for station keeping
                % dV_att = 30 * MAR_030 ; % [m/s] - Cumulative delta v for attitude control
                % 
                % orbiter.dV_disp = 250 * MAR_020 ; % [m/s] - Delta v for disposal of the orbiter
                % orbiter.dV_sk_o = dV_sk * 8 / 9 ; % When only orbiter ( 24 months )
                % orbiter.dV_att_o = dV_att * 8 / 9 ; % When only orbiter ( 24 months )
                % 
                % orbiter.dV_sk_cl = dV_sk / 9  ; % When clamped ( 3 months )
                % orbiter.dV_att_cl = dV_att / 9 ; % When clamped ( 192 months )
                % orbiter.dV_EOI = MAR_020 * 1100 + MAR_100 ; % [m/s] - Delta v for moon tour and Enceladus orbit inseriton (includes margin for moon approach navigation)
                % orbiter.dV_cap = 0.2477 * MAR_010 ;
                % orbiter.dV_int = 785.5 * MAR_010 + MAR_080 + MAR_090 ; % [m/s] - Delta v for interplanetary leg (includes margins on launcher dispersion and fly by)
                % orbiter.Isp_main = 321 ; % [s] - Specific impulse of main thrusters
                % orbiter.Isp_sk = 292 ; % [s] - Specific impulse of sk thrusters
                % orbiter.Isp_att = 292 ; % [s] - Specific impulse of attitude thrusters (same as sk)
                % 
                % orbiter.dry = 2359.1 * 1.2 ; % [kg] - Orbiter dry mass
                % lander.dry = 465.3 * 1.2 ; % [kg] - Lander dry mass
                % lander.hazard = 80 ; % [kg] - Hazard margin
                % 
                % % OVERALL PROPELLANT MASS - both modules
                % [ m_prop_orb, m_prop_lan ] = propellant_mass( orbiter, lander ) ;
                % 
                % fprintf( '-------------------------------- Old thrusters ------------------------------- \n' )
                % fprintf( 'Orbiter propellant mass: %f kg \n', m_prop_orb )
                % fprintf( 'Lander propellant mass: %f kg \n', m_prop_lan )
                % 
                % 
                % %% DATA - NEW CONFIGURATION (only orbiter)
                % 
                % % Lander
                % lander.dV_des = 244.6 * MAR_050 ; % [m/s] - Delta v for descending phase
                % lander.dV_att = 10 * MAR_030 ; % [m/s] - Delta v for attitude control during landing phase
                % lander.Isp_main = 228 ; % [s] - Specific impulse of main thrusters for landing
                % lander.Isp_att = 224 ; % [s] - Specific impulse of attitude thrusters for landing
                % 
                % % Orbiter
                % dV_sk = 370 * MAR_020 ; % [m/s] - Cumulative delta v for station keeping
                % dV_att = 30 * MAR_030 ; % [m/s] - Cumulative delta v for attitude control
                % 
                % orbiter.dV_disp = 250 * MAR_020 ; % [m/s] - Delta v for disposal of the orbiter
                % orbiter.dV_sk_o = dV_sk * 8 / 9 ; % When only orbiter ( 24 months )
                % orbiter.dV_att_o = dV_att * 8 / 9 ; % When only orbiter ( 24 months )
                % 
                % orbiter.dV_sk_cl = dV_sk / 9  ; % When clamped ( 3 months )
                % orbiter.dV_att_cl = dV_att / 9 ; % When clamped ( 192 months )
                % orbiter.dV_EOI = MAR_020 * 1100 + MAR_100 ; % [m/s] - Delta v for moon tour and Enceladus orbit inseriton (includes margin for moon approach navigation)
                % orbiter.dV_cap = 0.2477 * MAR_010 ;
                % orbiter.dV_int = 785.5 * MAR_010 + MAR_080 + MAR_090 ; % [m/s] - Delta v for interplanetary leg (includes margins on launcher dispersion and fly by)
                % orbiter.Isp_main = 329 ; % [s] - Specific impulse of main thrusters
                % orbiter.Isp_sk = 310 ; % [s] - Specific impulse of sk thrusters
                % orbiter.Isp_att = 227.5 ; % [s] - Specific impulse of attitude thrusters (same as sk)
                % 
                % orbiter.dry = 2359.1 * 1.2 ; % [kg] - Orbiter dry mass
                % lander.dry = 465.3 * 1.2 ; % [kg] - Lander dry mass
                % lander.hazard = 80 ; % [kg] - Hazard margin
                % 
                % % OVERALL PROPELLANT MASS - both modules
                % [ m_prop_orb, m_prop_lan ] = propellant_mass( orbiter, lander ) ;
                % 
                % fprintf( '---------------------------- New orbiter thrusters --------------------------- \n' )
                % fprintf( 'Orbiter propellant mass: %f kg \n', m_prop_orb )
                % fprintf( 'Lander propellant mass: %f kg \n', m_prop_lan )
                % 
                % 
                % %% DATA - NEW CONFIGURATION (only lander)
                % 
                % % Lander
                % lander.dV_des = 244.6 * MAR_050 ; % [m/s] - Delta v for descending phase
                % lander.dV_att = 10 * MAR_030 ; % [m/s] - Delta v for attitude control during landing phase
                % lander.Isp_main = 234 ; % [s] - Specific impulse of main thrusters for landing
                % lander.Isp_att = 227.5 ; % [s] - Specific impulse of attitude thrusters for landing
                % 
                % % Orbiter
                % dV_sk = 370 * MAR_020 ; % [m/s] - Cumulative delta v for station keeping
                % dV_att = 30 * MAR_030 ; % [m/s] - Cumulative delta v for attitude control
                % 
                % orbiter.dV_disp = 250 * MAR_020 ; % [m/s] - Delta v for disposal of the orbiter
                % orbiter.dV_sk_o = dV_sk * 8 / 9 ; % When only orbiter ( 24 months )
                % orbiter.dV_att_o = dV_att * 8 / 9 ; % When only orbiter ( 24 months )
                % 
                % orbiter.dV_sk_cl = dV_sk / 9  ; % When clamped ( 3 months )
                % orbiter.dV_att_cl = dV_att / 9 ; % When clamped ( 192 months )
                % orbiter.dV_EOI = MAR_020 * 1100 + MAR_100 ; % [m/s] - Delta v for moon tour and Enceladus orbit inseriton (includes margin for moon approach navigation)
                % orbiter.dV_cap = 0.2477 * MAR_010 ;
                % orbiter.dV_int = 785.5 * MAR_010 + MAR_080 + MAR_090 ; % [m/s] - Delta v for interplanetary leg (includes margins on launcher dispersion and fly by)
                % orbiter.Isp_main = 321 ; % [s] - Specific impulse of main thrusters
                % orbiter.Isp_sk = 292 ; % [s] - Specific impulse of sk thrusters
                % orbiter.Isp_att = 292 ; % [s] - Specific impulse of attitude thrusters (same as sk)
                % 
                % orbiter.dry = 2359.1 * 1.2 ; % [kg] - Orbiter dry mass
                % lander.dry = 465.3 * 1.2 ; % [kg] - Lander dry mass
                % lander.hazard = 80 ; % [kg] - Hazard margin
                % 
                % % OVERALL PROPELLANT MASS - both modules
                % [ m_prop_orb, m_prop_lan ] = propellant_mass( orbiter, lander ) ;
                % 
                % fprintf( '---------------------------- New lander thrusters ---------------------------- \n' )
                % fprintf( 'Orbiter propellant mass: %f kg \n', m_prop_orb )
                % fprintf( 'Lander propellant mass: %f kg \n', m_prop_lan )

%% DATA - NEW CONFIGURATION (both orbiter and lander)

% Lander
lander.dV_des = 244.6 * MAR_050 ; % [m/s] - Delta v for descending phase
lander.dV_att = 10 * MAR_030 ; % [m/s] - Delta v for attitude control during landing phase
lander.Isp_main = 234 ; % [s] - Specific impulse of main thrusters for landing
lander.Isp_att = 224 ; % [s] - Specific impulse of attitude thrusters for landing

% Orbiter
dV_sk = 2 * 212 * MAR_010 ; % [m/s] - Cumulative delta v for station keeping
%dV_att = 30 * MAR_030 ; % [m/s] - Cumulative delta v for attitude control

orbiter.dV_disp = 250 * MAR_020 ; % [m/s] - Delta v for disposal of the orbiter
orbiter.dV_sk_o = dV_sk * 8 / 9 ; % When only orbiter ( 24 months )
orbiter.m_att_science_o = ( 14.56 ) * MAR_attitude ; % When only orbiter ( 24 months )

%+ 54 * 189 / 216

orbiter.dV_sk_cl = dV_sk / 9  ; % When clamped ( 3 months )
orbiter.m_att_science_cl = ( 106.3 ) * MAR_attitude ; % When clamped ( 192 months ) + 20% margin
%+ 54 * 3 / 216 
orbiter.m_att_dsm_cl = 0.15 * MAR_attitude ; % When clamped ( 192 months ) + 20% margin
orbiter.m_att_moon_cl = 0.5 * MAR_attitude ; % During moon tour + 20% margin
orbiter.dV_EOI = ( 500 + 163 ) * MAR_010 + ( 111 + 370 + 142 + 53 + 105 ) * 1.25 + MAR_100 ; % [m/s] - Delta v for moon tour and Enceladus orbit inseriton (includes margin for moon approach navigation)
orbiter.dV_cap = 0.2477 * MAR_010 ;
orbiter.dV_int = 785.5 * MAR_010 + MAR_080 + MAR_090 ; % [m/s] - Delta v for interplanetary leg (includes margins on launcher dispersion and fly by)
orbiter.Isp_main = 329 ; % [s] - Specific impulse of main thrusters
orbiter.Isp_sk = 310 ; % [s] - Specific impulse of sk thrusters
orbiter.Isp_att = 224 ; % [s] - Specific impulse of attitude thrusters (same as sk)

orbiter.dry = 2390.3 * 1.2 ; % [kg] - Orbiter dry mass
lander.dry = 508.2 * 1.2 ; % [kg] - Lander dry mass
lander.hazard = 80 ; % [kg] - Hazard margin

% OVERALL PROPELLANT MASS - both modules
[ m_prop_orb, m_prop_lan, m_prop_main_orb, m_prop_sk_orb, m_prop_att_orb, m_prop_main_lan, m_prop_att_lan ] = propellant_mass( orbiter, lander ) ;

fprintf( '---------------------- New lander and orbiter thrusters ---------------------- \n' )
fprintf( 'Orbiter propellant mass: %f kg \n', m_prop_orb )
fprintf( 'Lander propellant mass: %f kg \n', m_prop_lan )


%% TANK SIZING - NEW CONFIGURATION - ORBITER
orbiter_data.prop_main = m_prop_main_orb ; % [kg] - Propellant mass required by main thrusters
orbiter_data.prop_sk = m_prop_sk_orb ; % [kg] - Propellant mass required by sk thrusters
orbiter_data.prop_att = m_prop_att_orb ; % [kg] - Propellant mass required by attitude thrusters

% Main thruster
main_thruster.Pc_mean = 11e5 ; % [ Pa ] - Mean combustion chamber pressure
main_thruster.of = 1.15 ; % [-] - Oxidizer to fuel ratio of main engine (mean value)

% SK thrusters
sk_thruster.Pc_mean = 16.5e5 ; % [ Pa ] - Mean combustion chamber pressure
sk_thruster.of = 0.85 ; % [-] - Oxidizer to fuel ratio of sk thrusters
sk_thruster.B = 4 ;% [ - ] - Blowdown ratio

% Attitude thruster
att_thruster.Pc_mean = 16.2e5 ; % [ Pa ] - Mean combustion chamber pressure
att_thruster.B = 4 ;% [ - ] - Blowdown ratio

% Tanks - Titanium
tank.rho = 2780 ; % [kg/m3] - Material density
tank.sigma = 950 ; % [MPa] - Material strenght
tank.temperature = 25 + 273.15 ; % [K] - Tank temperature (assumed to be at ambient temperature)

% Fuel - Hydrazine
rho_f = 1010 ; % [kg/m3] - Oxidizer density

% Pressurant - Helium
pressurant.R = 2077.3 ; % [J/kgK] - Specific gas constant
pressurant.gamma = 1.67 ; % Specific heat ratio of the gas
pressurant.Pi = 300e5 ; % [Pa] - Initial pressure of the pressurant gas

fprintf( '---------------------- Orbiter tank sizing ---------------------- \n' )
[ orb_tank_ox, orb_tank_fu, orb_tank_gas, V_gas_sec, V_prop_sec ] = multimode_architecture( main_thruster, sk_thruster, att_thruster, orbiter_data, tank, rho_f, pressurant )

%% TANK SIZING - NEW CONFIGURATION - LANDER
lander_data.prop_main = m_prop_main_lan ; % [kg] - Propellant mass required by main thrusters
lander_data.prop_att = m_prop_att_lan ; % [kg] - Propellant mass required by attitude thrusters

% Main thruster
main_thruster.Pc_mean = 5.5e5 ; % [ Pa ] - Mean combustion chamber pressure
main_thruster.B = 4 ; % [-] - Blowdown ratio of main engine (mean value)

% Attitude thruster
att_thruster.Pc_mean = 16.2e5 ; % [ Pa ] - Mean combustion chamber pressure
att_thruster.B = 4 ; % [ - ] - Blowdown ratio of attitude thrusters

% Tanks - Titanium
tank.rho = 2780 ; % [kg/m3] - Material density
tank.sigma = 950 ; % [MPa] - Material strenght
tank.temperature = 25 + 273.15 ; % [K] - Tank temperature (assumed to be at ambient temperature)

% Fuel - Hydrazine
rho_f = 1010 ; % [kg/m3] - Oxidizer density

% Pressurant - Helium
pressurant.R = 2077.3 ; % [J/kgK] - Specific gas constant
pressurant.gamma = 1.67 ; % Specific heat ratio of the gas
pressurant.Pi = 300e5 ; % [Pa] - Initial pressure of the pressurant gas

fprintf( '---------------------- Lander tank sizing ---------------------- \n' )
[ lan_tank_prop, lan_tank_gas, V_prop_att, V_gas_att ] = dualmode_lander( main_thruster, att_thruster, lander_data, tank, rho_f, pressurant )

%% FIRING TIME

% ORBITER
mdot_main_orb = 180e-3 ; % [kg/s] - Maximum flow rate
t_firing_main_orb = m_prop_main_orb / mdot_main_orb  % [s] - Total firing time of primary propulsion

% LANDER
mdot_main_lan = 98e-3 ; % [kg/s] - Maximum flow rate
t_firing_main_lan = m_prop_main_lan / mdot_main_lan  % [s] - Total firing time of primary propulsion
