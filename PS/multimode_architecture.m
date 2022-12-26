function [ tank_ox, tank_fu, tank_gas ] = multimode_architecture( main_thruster, sk_thruster, att_thruster, orbiter_data, tank, rho_f, pressurant )

% main_thruster, sk_thruster, att_thruster structures:
% .Pc_mean  [Pa] - Mean combustion chamber pressure
% .of  Oxidizer to fuel ratio (not for the attitude thrusters)
% .B  Blowdown ratio (not needed for main thrusters)

% orbiter_data structure:
% .prop_main [kg] - Propellant mass required by main thrusters (already with 2.5% margin)
% .prop_sk  [kg] - Propellant mass required by sk thrusters (already with 2.5% margin)
% .prop_att  [kg] - Propellant mass required by sttitude thrusters (already with 20% margin)

% tank structure
% .rho  [kg/m3] - Density of material
% .sigma [MPa] - Material strenght

% rho_f  [kg/m3] - Fuel density

% pressurant structure
% .R  [J/kgK] - Specific gas constant
% .gamma  [-] -Specific heat ratio
% .Pi  [Pa] - Initial pressure of the pressurant gas

%% DATA

% General data
m_prop_main = orbiter_data.prop_main ; % [kg] - Propellant mass required by main thrusters
m_prop_sk = orbiter_data.prop_sk ; % [kg] - Propellant mass required by sk thrusters
m_prop_att = orbiter_data.prop_att ; % [kg] - Propellant mass required by attitude thrusters

dP_feed = 40e3 ; % [ Pa ] - Feeding system losses

% Main thruster
Pc_main = main_thruster.Pc_mean ; % [ Pa ] - Mean combustion chamber pressure
dP_inj_main = 0.3 * Pc_main ; % [ Pa ] - Injection losses
OF_main = main_thruster.of ; % [-] - Oxidizer to fuel ratio of main engine

% SK thrusters
Pc_sk = sk_thruster.Pc_mean ; % [ Pa ] - Mean combustion chamber pressure
dP_inj_sk = 0.3 * Pc_sk; % [ Pa ] - Injection losses
OF_sk = sk_thruster.of ; % [-] - Oxidizer to fuel ratio of sk thrusters
B_sk = sk_thruster.B ;% [ - ] - Blowdown ratio

% Attitude thruster
Pc_att = att_thruster.Pc_mean ; % [ Pa ] - Mean combustion chamber pressure
dP_inj_att = 0.3 * Pc_att ; % [ Pa ] - Injection losses
B_att = att_thruster.B ;% [ - ] - Blowdown ratio

% Tanks
rho_tank = tank.rho ; % [kg/m3] - Material density
sigma = tank.sigma * 1e6 ; % [Pa] - Material strenght
T_tank = tank.temperature ; % [K] - Tank temperature

% Pressurant
R_gas = pressurant.R ; % [J/kgK] - Specific gas constant
gamma = pressurant.gamma ; % Specific heat ratio of the gas
Pi_gas = pressurant.Pi ; % [Pa] - Initial pressure of the pressurant gas

%% Preliminary computations
m_fu_main = m_prop_main / ( 1 + OF_main ) * 1.01 ; % [kg] - Fuel mass + ullage margin (1%)
m_ox_main = ( m_prop_main - m_fu_main ) * 1.01 ; % [kg] - Oxidizer mass + ullage margin (1%)
V_fu_main = m_fu_main / rho_f ; % [m3] - Volume of the fuel
V_ox_main = V_fu_main ; % [m3] - Volume of the oxidizer (imposed to be equal)

m_fu_sk = m_prop_sk / ( 1 + OF_sk ) * 1.01 ; % [kg] - Fuel mass + ullage margin (1%)
m_ox_sk = ( m_prop_sk - m_fu_sk ) * 1.01 ; % [kg] - Oxidizer mass + ullage margin (1%)
V_fu_sk = m_fu_sk / rho_f ; % [m3] - Volume of the fuel
V_ox_sk = V_fu_sk ; % [m3] - Volume of the oxidizer (imposed to be equal)

V_prop_att = m_prop_att / rho_f ; % [m3] - Volume of the propellant 

Pi_prop_main = Pc_main + dP_inj_main + dP_feed ; % [ Pa ] - Initial propellant pressure (thrusters)

Pi_prop_sk = Pc_sk + dP_inj_sk + dP_feed ; % [ Pa ] - Initial propellant pressure (sk thrusters)

Pi_prop_att = Pc_att + dP_inj_main + dP_inj_att ; % [ Pa ] - Initial propellant pressure 
Pi_prop_tank = max( [ Pi_prop_main, Pi_prop_sk, Pi_prop_att ] ) ; % [Pa] - Initial propellant (both ox and fu) tank pressure (maximum between the required one) - it has to be the same for all thrusters as we want to have only one common tank

%% PRESSURE REGULATED - BIPROPELLANT THRUSTER (MAIN)

% Oxidizer
m_gas_ox_main = ( Pi_prop_tank * V_ox_main ) / ( R_gas * T_tank ) * gamma / ...  % [ kg ] - Pressurant mass for the oxidizer
           ( 1 - Pi_prop_tank / Pi_gas ) ;

V_gas_ox_main = ( m_gas_ox_main * R_gas * T_tank ) / Pi_gas ;  % [ m3 ] - Pressurant volume for the oxidizer (coincides with initial pressurant tank volume)

% Fuel
m_gas_fu_main = ( Pi_prop_tank * V_fu_main ) / ( R_gas * T_tank ) * gamma / ...  % [ kg ] - Pressurant mass for the orbiter
           ( 1 - Pi_prop_tank / Pi_gas ) ;

V_gas_fu_main = ( m_gas_fu_main * R_gas * T_tank ) / Pi_gas ;  % [ m3 ] - Pressurant volume needed for the fuel (coincides with initial pressurant tank volume)

m_gas_main = m_gas_fu_main + m_gas_ox_main ; % [kg] - Pressurant mass needed for main propellant
V_gas_main = V_gas_fu_main + V_gas_ox_main ; % [m3] - Pressurant volume needed for main propellant

%% BLOWDOWN - BIPROPELLANT THRUSTER (SK)
V_gas_sk = ( V_ox_sk + V_fu_sk ) / ( B_sk - 1 ) ; % [m3] - Initial volume of pressurant gas
m_gas_sk = ( Pi_gas * V_gas_sk ) / ( R_gas * T_tank ) ; % [kg] - Pressurant mass 

%% BLOWDOWN - MONOPROPELLANT THRUSTERS (ATTITUDE)
V_gas_att = V_prop_att / ( B_att - 1 ) ; % [m3] - Initial volume of pressurant gas
m_gas_att = ( Pi_gas * V_gas_att ) / ( R_gas * T_tank ) ; % [kg] - Pressurant mass 

%% CUMULATIVE MASSES
m_gas = ( m_gas_main + m_gas_sk + m_gas_att ) * 1.2 ; % [kg] - Overall pressurant mass
m_ox = m_ox_main + m_ox_sk ; % [kg] - Overall oxidizer mass
m_fu = m_fu_main + m_fu_sk + m_prop_att ; % [kg] - Overall fuel mass

V_gas = V_gas_main + V_gas_sk + V_gas_att ; % [m3] - Overall pressurant mass 
V_ox = ( V_ox_main + V_ox_sk ) * 1.1 ; % [m3] - Overall oxidizer volume + 10% of unsable volume margin
V_fu = ( V_fu_main + V_fu_sk + V_prop_att ) * 1.1 ; % [m3] - Overall fuel volume + 10% of unsable volume margin

%% TANK SIZING - Assuming that oxidizer, fuel and pressurant are in three separate tanks

% Pressurant tank - spherical
r_tank_gas = ( 3 * V_gas / ( 4 * pi ) ) ^ ( 1 / 3 ) ; % [m] - Radius of the pressurant tank
t_tank_gas = Pi_gas * r_tank_gas / ( 2 * sigma ) ; % [m] - Pressurant tank thickness
m_tank_gas = 4 / 3 * rho_tank * pi * ( ( r_tank_gas + t_tank_gas ) ^ 3 - r_tank_gas ^ 3  ) ; % [kg] - Pressurant tank mass

% Oxidizer tank - cylindrical with hemispherical caps
r_tank_ox = ( 3 * V_ox / ( 7 * pi ) ) ^ ( 1 / 3 ) ; % [m] - Radius of the oxidizer tank
h_tank_ox = r_tank_ox ; % [m] - Height of the oxidizer tank (assumed to be equal to the radius as it is the best compromise for structural reason)
t_tank_ox = Pi_prop_tank * r_tank_ox / sigma ; % [m] - Oxidizer tank thickness
m_tank_ox = 4 / 3 * rho_tank * pi * ( ( r_tank_ox + t_tank_ox ) ^ 3 - r_tank_ox ^ 3  ) ; % [kg] - Oxidizer tank mass

% Fuel tank - cylindrical with hemispherical caps
r_tank_fu = ( 3 * V_fu / ( 7 * pi ) ) ^ ( 1 / 3 ) ; % [m] - Radius of the oxidizer tank
h_tank_fu = r_tank_fu ; % [m] - Height of the oxidizer tank (assumed to be equal to the radius as it is the best compromise for structural reason)
t_tank_fu = Pi_prop_tank * r_tank_fu / sigma ; % [m] - Oxidizer tank thickness
m_tank_fu = 4 / 3 * rho_tank * pi * ( ( r_tank_fu + t_tank_fu ) ^ 3 - r_tank_fu ^ 3  ) ; % [kg] - Oxidizer tank mass

%% Define sizing outputs

tank_ox.ox_mass = m_ox ;
tank_ox.ox_volume = V_ox ;
tank_ox.ox_tank_radius = r_tank_ox ;
tank_ox.ox_tank_height = h_tank_ox ;
tank_ox.ox_tank_thickness = t_tank_ox * 1e3 ; % [mm]
tank_ox.ox_tank_mass = m_tank_ox ;

tank_fu.fu_mass = m_fu ;
tank_fu.fu_volume = V_fu ;
tank_fu.fu_tank_radius = r_tank_fu ;
tank_fu.fu_tank_height = h_tank_fu ;
tank_fu.fu_tank_thickness = t_tank_fu * 1e3 ; % [mm]
tank_fu.fu_tank_mass = m_tank_fu ;

tank_gas.gas_mass = m_gas ;
tank_gas.gas_volume = V_gas ;
tank_gas.gas_tank_radius = r_tank_gas ;
tank_gas.gas_tank_thickness = t_tank_gas * 1e3 ; % [mm]
tank_gas.gas_tank_mass = m_tank_gas ;
