function [ tank_prop, tank_gas ] = dualmode_lander( main_thruster, att_thruster, lander_data, tank, rho_f, pressurant )

% main_thruster, sk_thruster, att_thruster structures:
% .Pc_mean  [Pa] - Mean combustion chamber pressure
% .B  Blowdown ratio (not needed for main thrusters)

% lander_data structure:
% .prop_main [kg] - Propellant mass required by main thrusters (already with 2.5% margin)
% .prop_att  [kg] - Propellant mass required by sttitude thrusters (already with 2.5% margin)

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
m_prop_main = lander_data.prop_main ; % [kg] - Propellant mass required by main thrusters
m_prop_att = lander_data.prop_att ; % [kg] - Propellant mass required by attitude thrusters

dP_feed = 40e3 ; % [ Pa ] - Feeding system losses

% Main thruster
Pc_main = main_thruster.Pc_mean ; % [ Pa ] - Mean combustion chamber pressure
dP_inj_main = 0.3 * Pc_main ; % [ Pa ] - Injection losses
B_main = main_thruster.B ;% [ - ] - Blowdown ratio

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
m_prop_main = m_prop_main * 1.01 ; % [kg] - Oxidizer mass + ullage margin (1%)
V_prop_main = m_prop_main / rho_f ; % [m3] - Volume of the propellant for main thrusters

m_prop_att = m_prop_att * 1.01 ; % [kg] - Oxidizer mass + ullage margin (1%)
V_prop_att = m_prop_att / rho_f ; % [m3] - Volume of the propellant for attitude 

m_prop = m_prop_att + m_prop_main ; % [kg] - Propellant mass
V_prop = V_prop_att + V_prop_main ; % [m3] - Volume of the propellant 

Pi_prop_main = Pc_main + dP_inj_main + dP_feed ; % [ Pa ] - Initial propellant pressure (thrusters)
Pi_prop_att = Pc_att + dP_inj_main + dP_inj_att ; % [ Pa ] - Initial propellant pressure 

Pi_prop_tank = max( [ Pi_prop_main, Pi_prop_att ] ) ; % [Pa] - Initial propellant (both ox and fu) tank pressure (maximum between the required one) - it has to be the same for all thrusters as we want to have only one common tank

%% BLOWDOWN - MAIN MONOPROPELLANT THRUSTERS
V_gas_main = V_prop_main / ( B_main - 1 ) ; % [m3] - Initial volume of pressurant gas
m_gas_main = ( Pi_gas * V_gas_main ) / ( R_gas * T_tank ) ; % [kg] - Pressurant mass 

%% BLOWDOWN - MONOPROPELLANT THRUSTERS (ATTITUDE)
V_gas_att = V_prop_att / ( B_att - 1 ) ; % [m3] - Initial volume of pressurant gas
m_gas_att = ( Pi_gas * V_gas_att ) / ( R_gas * T_tank ) ; % [kg] - Pressurant mass 

%% CUMULATIVE MASSES
m_gas = ( m_gas_main + m_gas_att ) * 1.2 ; % [kg] - Overall pressurant mass
V_gas = V_gas_main + V_gas_att ; % [m3] - Overall pressurant mass 

%% TANK SIZING - Assuming that oxidizer, fuel and pressurant are in three separate tanks

% Pressurant tank - spherical
r_tank_gas = ( 3 * V_gas / ( 4 * pi ) ) ^ ( 1 / 3 ) ; % [m] - Radius of the pressurant tank
t_tank_gas = Pi_gas * r_tank_gas / ( 2 * sigma ) ; % [m] - Pressurant tank thickness
m_tank_gas = 4 / 3 * rho_tank * pi * ( ( r_tank_gas + t_tank_gas ) ^ 3 - r_tank_gas ^ 3  ) ; % [kg] - Pressurant tank mass

% Fuel tank - cylindrical with hemispherical caps
r_tank_prop = ( 3 * V_prop / ( 7 * pi ) ) ^ ( 1 / 3 ) ; % [m] - Radius of the oxidizer tank
h_tank_prop = r_tank_prop ; % [m] - Height of the oxidizer tank (assumed to be equal to the radius as it is the best compromise for structural reason)
t_tank_prop = Pi_prop_tank * r_tank_prop / sigma ; % [m] - Oxidizer tank thickness
m_tank_prop = 4 / 3 * rho_tank * pi * ( ( r_tank_prop + t_tank_prop ) ^ 3 - r_tank_prop ^ 3  ) ; % [kg] - Oxidizer tank mass

%% Define sizing outputs

tank_prop.prop_mass = m_prop ;
tank_prop.prop_volume = V_prop ;
tank_prop.prop_tank_radius = r_tank_prop ;
tank_prop.prop_tank_height = h_tank_prop ;
tank_prop.prop_tank_thickness = t_tank_prop * 1e3 ; % [mm]
tank_prop.prop_tank_mass = m_tank_prop ;

tank_gas.gas_mass = m_gas ;
tank_gas.gas_volume = V_gas ;
tank_gas.gas_tank_radius = r_tank_gas ;
tank_gas.gas_tank_thickness = t_tank_gas * 1e3 ; % [mm]
tank_gas.gas_tank_mass = m_tank_gas ;