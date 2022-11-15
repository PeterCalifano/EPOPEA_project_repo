function [ M_PS, sizing_ox, sizing_fu, sizing_gas ] = bipropellant_DM( oxidizer, fuel, pressurant, tank, thruster )
% BIPROPELLANT sizing of tanks in dual mode
% A part of mp is feed with constan pressure (mp_reg) -> feeding MAIN ENGINE
% The rest of mp is fed in blowdown (mp_blow) -> feeding SMALLER THRUSTERS

% INPUTS
%
% oxidizer      [ struc ] - Structure with all propellant properties:
%                               oxidizer.rho : density of propellant                 [ kg/m^3 ]          
%                               oxidizer.OFratio_reg : O/F ratio of Main Engines     [ - ]
%                               oxidizer.OFratio_blow : O/F ratio of Small Thrusters [ - ]
%                               oxidizer.prop_mass_reg : mass of propellant to feed with regulated pressure [ kg ]
%                               oxidizer.prop_mass_blow : mass of propellant to feed with blow-down mode [ kg ]
%                               oxidizer.type : geometry of oxidizer tank 
%                                               ('Spherical' or 'Cylindrical')      
%
% fuel          [ struc ] - Structure with all propellant properties:
%                               fuel.rho : density of propellant                  [ kg/m^3 ]   

% pressurant    [ struc ] - Structure with all pressurant properties:
%                               pressurant.feed_system.B : Blowdown ratio (only   [ - ]
%                                                          if blowdown system)
%                               pressurant.gamma : specific heat ratio            [ - ]
%                               pressurant.R : specific gas constant              [ J/kgK ]
%                               pressurant.pressure : initial pressure of         [ Pa ]
%                                                     pressurant (consider the one for the initial pressure regulated mode)  

% tank          [ struc ] - Structure with all tank properties
%                               tank.rho : tank material density                  [ kg/m^3 ]
%                               tank.sigma : tank material strenght               [ Pa ]
%                               tank.ratio : cylindrical intersection radiut-to-height ratio    [ m ]
%                               tank.temperature : tank temperature (if not given [ K ]
%                                                  will be considered as 10 deg)               
                        %                               tank.type : tank geometry                         [ str ]
                        %                                       - 'Spherical' 
                        %                                       - 'Cylindrical' -> cylindrical with
                        %                                                          spherical caps
%
% thruster      [ struc ] - Structure with all tank properties
%                               thruster.mass_ME : Main Engine mass                     [ kg ]
%                               thruster.number_ME : number of Main Engines             [ - ]
%                               thruster.mass_ST : Small Thrusters mass                 [ kg ]
%                               thruster.number_ST : number of Small Thrusters          [ - ]
%                               thruster.chamber_pressure_reg : combustion chamber pressure in ME           [ Pa ]
%                               thruster.chamber_pressure_blow : combustion chamber pressure in ST          [ Pa ]

% OUTPUTS
%
% M_PS          [ 1 ]     - Total mass of the propulsion subsystem                [ kg ]
%
% sizing_ox     [ struc ] - Structure with tank sizing:
%                               sizing_ox.single : tank properties (when single  [ struc ] 
%                                                   tank for ox and pressurant is 
%                                                   used)         
%                                   - .type                                       [ str ]
%                                   - .mass_press                                 [ kg ]
%                                   - .volume_press                               [ m^3 ]
%                                   - .mass_ox                                    [ kg ]
%                                   - .volume_ox                                  [ m^3 ]
%                                   - .radius                                     [ m ]
%                                   - .thickness                                  [ m ]
%                                   - .tank_mass                                  [ kg ]
%                                   - .tank_volume                                [ m^3 ]
%
%                               sizing_ox.type                                    [ str ]
%                               sizing_ox.mass                                    [ kg ]
%                               sizing_ox.mass_press                              [ kg ]
%                               sizing_ox.volume                                  [ m^3 ]
%                               sizing_ox.volume_press                            [ m^3 ]
%                               sizing_ox.radius                                  [ m ]
%                               sizing_ox.radius_press                            [ m ]
%                               sizing_ox.thickness                               [ m ]
%                               sizing_ox.thickness_press                         [ m ]
%                               sizing_ox.tank_mass                               [ kg ]
%                               sizing_ox.tank_mass_press                         [ kg ]

% sizing_fu     [ struc ] - Structure with fuel tank sizing:
%                               sizing_fu.single : tank properties (when single  [ struc ] 
%                                                   tank for ox and pressurant is 
%                                                   used)         
%                                   - .type                                       [ str ]
%                                   - .mass_press                                 [ kg ]
%                                   - .volume_press                               [ m^3 ]
%                                   - .mass_fu                                    [ kg ]
%                                   - .volume_fu                                  [ m^3 ]
%                                   - .radius                                     [ m ]
%                                   - .thickness                                  [ m ]
%                                   - .tank_mass                                  [ kg ]
%                                   - .tank_volume                                [ m^3 ]
%
%                               sizing_fu.type                                    [ str ]
%                               sizing_fu.mass                                    [ kg ]
%                               sizing_fu.mass_press                              [ kg ]
%                               sizing_fu.volume                                  [ m^3 ]
%                               sizing_fu.volume_press                            [ m^3 ]
%                               sizing_fu.radius                                  [ m ]
%                               sizing_fu.radius_press                            [ m ]
%                               sizing_fu.thickness                               [ m ]
%                               sizing_fu.thickness_press                         [ m ]
%                               sizing_fu.tank_mass                               [ kg ]
%                               sizing_fu.tank_mass_press                         [ kg ]

% sizing_gas    [ struc ] - Structure with pressurant tank sizing:
%                               sizing_gas.type                                   [ str ]
%                               sizing_gas.mass                                   [ kg ]
%                               sizing_gas.volume                                 [ m^3 ]
%                               sizing_gas.radius                                 [ m ]
%                               sizing_gas.thickness                              [ m ]
%                               sizing_gas.tank_mass                              [ kg ]
%                           NOTE: This structure is given as output only in
%                           case of n=3 tanks.

% ---

%% DATA
mp_reg = oxidizer.prop_mass_reg ;                                                 % [ kg ] - Propellant mass for pressure regulated feeding 
mp_reg = mp_reg * 1.12 ;                                                          % [ kg ] Margin: MAR-CP-010 + 3% of unsable mass
mp_blow = oxidizer.prop_mass_blow ;                                               % [ kg ] - Propellant mass for blow-down feeding 
mp_blow = mp_blow * 1.12 ;                                                        % [ kg ] Margin: MAR-CP-010 + 3% of unsable mass
  
OF_reg = oxidizer.OFratio_reg ;                                                   % [ - ] - Oxidizer-to-fuel ratio of Main Engine
OF_blow = oxidizer.OFratio_blow ;                                                 % [ - ] - Oxidizer-to-fuel ratio of smaller thrusters

mf_reg = mp_reg / ( 1 + OF_reg ) ;                                                    % [ kg ] - Fuel mass (pressure regulated portion)
mf_reg = 1.01 * mf_reg ;                                                              % [ kg ] - Margin: ullage
mf_blow = mp_blow / ( 1 + OF_blow ) ;                                                 % [ kg ] - Fuel mass (blow-down portion)
mf_blow = 1.01 * mf_blow ;                                                            % [ kg ] - Margin: ullage

mox_reg = ( OF_reg * mp_reg ) / ( OF_reg + 1 ) ;                                      % [ kg ] - Oxidizer mass (pressure regulated portion)
mox_reg = 1.01 * mox_reg ;                                                            % [ kg ] - Margin: ullage
mox_blow = ( OF_blow * mp_blow ) / ( OF_blow + 1 ) ;                                  % [ kg ] - Oxidizer mass (blow-down portion)
mox_blow = 1.01 * mox_blow ;                                                          % [ kg ] - Margin: ullage

rho_ox = oxidizer.rho ;                                                           % [ kg/m^3 ] - Oxidizer density
rho_fu = fuel.rho ;                                                               % [ kg/m^3 ] - Fuel density

rho_p_reg = mp_reg / ( mf_reg / rho_fu + mox_reg / rho_ox ) ;                     % [ kg/m^3 ] - Mean propellant density (pressure regulated portion)
rho_p_blow = mp_blow / ( mf_blow / rho_fu + mox_blow / rho_ox ) ;                 % [ kg/m^3 ] - Mean propellant density (blowdown portion)

Vox_reg = mox_reg / rho_ox ;                                                        % [ m^3 ] - Oxidizer volume (pressure regulated portion)
Vf_reg = mf_reg / rho_fu ;                                                          % [ m^3 ] - Fuel volume (pressure regulated portion)
Vox_blow = mox_blow / rho_ox ;                                                        % [ m^3 ] - Oxidizer volume (blow-down portion)
Vf_blow = mf_blow / rho_fu ;                                                          % [ m^3 ] - Fuel volume (blow-down portion)

Vp_reg = mp_reg / rho_p_reg ;                                                      % [ m^3 ] - Propellant volume (pressure regulated portion)
Vp_blow = mp_blow / rho_p_blow ;                                                   % [ m^3 ] - Propellant volume (blow-down portion)

% If the initial pressure of the propellant is not given as input, it is
% retrieved from pressure losses and desired combustion chamber pressure 
if isfield( thruster, 'chamber_pressure' )
    Pc_reg = thruster.chamber_pressure_reg ;                                      % [ Pa ] - Combustion chamber pressure of the Main Engine
    Pc_blow = thruster.chamber_pressure_blow ;                                     % [ Pa ] - Combustion chamber pressure of the Small Thrusters
else
    Pc_reg = 20e5 ;      
    Pc_blow = 10e5 ;                                                              % (The lower N thrusters usually require a lower Pc)   
end

dP_inj_reg = 0.3 * Pc_reg ;                                                              % [ Pa ] - Injection losses
dP_feed = 40e3 ;                                                                         % [ Pa ] - Feeding system losses
Pi_p_reg = Pc_reg + dP_inj_reg + dP_feed ;                                               % [ Pa ] - Propellant tank initial pressure
dP_inj_blow = 0.3 * Pc_blow ;                                                            % [ Pa ] - Injection losses                                                                       % [ Pa ] - Feeding system losses
Pi_p_blow = Pc_blow + dP_inj_blow + dP_feed ;                                            % [ Pa ] - Propellant tank initial pressure
Pi_prop = max(Pi_p_blow,Pi_p_reg);                                                       % Take maximum pressure to pressurize the propellant tanks

% Pressurant is common for both modes
gamma = pressurant.gamma ;                                                        % [ - ] - Specific heat ratio of pressurant
R = pressurant.R ;                                                                % [ J/kgK ] - Specific gas constant fo pressurant
Pi_gas = pressurant.pressure ;                                                    % [ Pa ] - Initial pressurant pressure

rho_tank = tank.rho ;                                                             % [ kg/m^3] - Tank material density
sigma_tank = tank.sigma ;                                                         % [ Pa ] - Tank material strength

% If the tank temperature is not given as input, the ambient temperature 
% (10 degrees) is considered
if isfield( tank, 'temperature' )
    T_tank = tank.temperature ;                                                   % [ k ] - Tank temperature
else
    T_tank = 293 ;                                                               
end
               

%% PRESSURE REGULATED
% Oxidizer
m_gas_ox = ( Pi_p_reg * Vox_reg ) / ( R * T_tank ) * gamma / ...           % [ kg ] - Pressurant mass
        ( 1 - Pi_p_reg / Pi_gas ) ;                             
m_gas_ox = m_gas_ox * 1.2 ;                                                % [ kg ] - Margin: MAR-MAS-090

V_gas_ox_reg = ( m_gas_ox * R * T_tank ) / Pi_gas ;                        % [ m^3 ] - Pressurant volume (coincides with initial pressurant tank volume)
V_gas_ox_reg = 1.01 * V_gas_ox_reg ;                                       % [ m^3 ] - Margin: bladder 

% Fuel
m_gas_fu = ( Pi_p_reg * Vf_reg ) / ( R * T_tank ) * gamma / ...            % [ kg ] - Pressurant mass
        ( 1 - Pi_p_reg / Pi_gas ) ;                             
m_gas_fu = m_gas_fu * 1.2 ;                                                % [ kg ] - Margin: MAR-MAS-090      

V_gas_fu_reg = ( m_gas_fu * R * T_tank ) / Pi_gas ;                        % [ m^3 ] - Pressurant volume (coincides with initial pressurant tank volume)
V_gas_fu_reg = 1.01 * V_gas_fu_reg ;                                       % [ m^3 ] - Margin: bladder 

m_gas_reg = m_gas_ox + m_gas_fu ;                                          % Total pressurant mass and volume
                                                                           % required for the regulated portion
V_gas_reg = V_gas_ox_reg + V_gas_fu_reg ;

%% BLOWDOWN
B = pressurant.B ;                                                         % [ - ] - Blowdown ratio

% Check on B value
P_c_min = min(Pc_reg,Pc_blow) ;
Pf_gas = P_c_min + dP_feed ;
B_max = Pi_gas / Pf_gas ;
if B > B_max
    error( 'Blowdown ratio exceeds the maximum value: %d \n', B_max )
end

% Oxidizer
V_gas_ox_blow = Vox_blow / ( B - 1) ;                                      % [ m^3 ] - Pressurant volume for oxidizer
V_gas_ox_blow = V_gas_ox_blow * 1.01 ;                                     % [ m^3 ] - Margin: baldder volume

m_gas_ox = ( Pi_gas * V_gas_ox_blow ) / ( R * T_tank ) ;                   % [ kg ] - Pressurant mass for oxidizer
m_gas_ox = m_gas_ox * 1.2 ;                                                % [ kg ] - Margin: MAR-MAS-090

% Fuel
V_gas_fu_blow = Vf_blow / ( B - 1) ;                                       % [ m^3 ] - Pressurant volume for fuel
V_gas_fu_blow = V_gas_fu_blow * 1.01 ;                                     % [ m^3 ] - Margin: baldder volume

m_gas_fu = ( Pi_gas * V_gas_fu_blow ) / ( R * T_tank ) ;                   % [ kg ] - Pressurant mass for fuel
m_gas_fu = m_gas_fu * 1.2 ;                                                % [ kg ] - Margin: MAR-MAS-090

m_gas_blow = m_gas_ox + m_gas_fu ;
V_gas_blow = V_gas_ox_blow + V_gas_fu_blow ;

%% UNIFY PRESSURE-REGULATED AND BLOWDOWN
% The same mass and tank of pressurant and propellant will be used, even if
% switching from one mode to the other
m_gas = m_gas_blow + m_gas_reg;
V_gas = V_gas_blow + V_gas_reg;
V_gas_ox = V_gas_ox_blow + V_gas_ox_reg;
V_gas_fu = V_gas_fu_blow + V_gas_fu_reg;
m_ox = mox_blow + mox_reg;
m_fu = mf_blow + mf_reg;
V_ox = Vox_blow + Vox_reg;
V_fu = Vf_blow + Vf_reg;

%% TANK SIZING
% 3 tanks: 1 for NTO, 1 for Hydrazine, 1 for pressurant

switch oxidizer.type
case 'Spherical'
    % Oxidizer tank
    r_tank_ox = ( ( 3 * V_ox ) / ( 4 * pi ) ) ^ ( 1 / 3 );                 % [ m ] - Propellant tank radius
    t_tank_ox = ( Pi_prop * r_tank_ox ) / ( 2 * sigma_tank ) ;             % [ m ] - Propellant tank thickness
    m_tank_ox = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_ox + ...            % [ kg ] - Propellant tank mass
                  t_tank_ox) ^ 3 - r_tank_ox ^ 3 ) ;

    % Fuel tank
    r_tank_fu = ( ( 3 * V_fu ) / ( 4 * pi ) ) ^ ( 1 / 3 );                 % [ m ] - Propellant tank radius
    t_tank_fu = ( Pi_prop * r_tank_fu ) / ( 2 * sigma_tank ) ;             % [ m ] - Propellant tank thickness
    m_tank_fu = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_fu + ...            % [ kg ] - Propellant tank mass
                  t_tank_fu) ^ 3 - r_tank_fu ^ 3 ) ;

case 'Cylindrical'
    % Oxidizer tank
    ratio = tank.ratio ;                                                   % [ m ] - Radius-to-height ratio of cylindrical interection - NOTE: we assumed both propellant and pressurant tanks to have same height
    % height = radius / ratio

    r_function_ox = @( x ) V_ox - ( 4 * pi / 3 ) * x ^ 3 - pi * ...       % Function used to compute the radius of the tank knowing the height
                             x / ratio * x ^ 2 ; 
    r_guess_ox = 0.3 * ratio ;
    r_tank_ox = fzero( r_function_ox, r_guess_ox ) ;   

    h_tank = r_tank_ox / ratio ;

    sizing_ox.heigth = h_tank ;                                    % Save in output structure

    t_tank_ox = ( Pi_prop * r_tank_ox ) / sigma_tank ;                    % [ m ] - Tank thickness
   
    m_tank_ox = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_ox + ...         % [ kg ] - Tank mass
                 t_tank_ox ) ^ 3 - r_tank_ox ^ 3 ) +  pi * h_tank * ...
                 ( ( r_tank_ox + t_tank_ox ) ^ 2 - r_tank_ox ^ 2 ) ) ;

    % Fuel tank
    r_function_fu = @( x ) V_fu - ( 4 * pi / 3 ) * x ^ 3 - pi * ...       % Function used to compute the radius of the tank knowing the height
                             x / ratio * x ^ 2 ; 
    r_guess_fu = 0.3 * ratio ;
    r_tank_fu = fzero( r_function_fu, r_guess_fu ) ;   

    h_tank = r_tank_fu / ratio ;
    sizing_fu.heigth = h_tank ;                                    % Save in output structure

    t_tank_fu = ( Pi_prop * r_tank_fu ) / sigma_tank ;                    % [ m ] - Tank thickness
   
    m_tank_fu = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_fu + ...         % [ kg ] - Tank mass
                 t_tank_fu ) ^ 3 - r_tank_fu ^ 3 ) +  pi * h_tank * ...
                 ( ( r_tank_fu + t_tank_fu ) ^ 2 - r_tank_fu ^ 2 ) ) ;

end

% For pressurant only a spherical tank is considered
r_tank_gas = ( ( 3 * ( V_gas_fu + V_gas_ox ) ) / ( 4 * pi ) ) ^ ...   % [ m ] - Pressurant tank radius
    ( 1 / 3 );
t_tank_gas = ( Pi_gas * r_tank_gas ) / ( 2 * sigma_tank ) ;           % [ m ] - Pressurant anks thickness
m_tank_gas = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_gas + ...         % [ kg ] - Pressurant tank mass
    t_tank_gas ) ^ 3 - r_tank_gas ^ 3 ) ;


m_tank = m_tank_gas + m_tank_ox + m_tank_fu ;                                 % [ kg ] - Total tank mass

%% RESULTS
sizing_ox.type = oxidizer.type ;                                              % Save in output structure
sizing_ox.mass = m_ox ;                                                       % Save in output structure
sizing_ox.volume = V_ox ;                                                     % Save in output structure
sizing_ox.radius = r_tank_ox ;                                                % Save in output structure
sizing_ox.thickness = t_tank_ox ;                                             % Save in output structure
sizing_ox.tank_mass = m_tank_ox ;                                             % Save in output structure

sizing_fu.type = oxidizer.type ;                                                  % Save in output structure
sizing_fu.mass = m_fu ;                                                       % Save in output structure
sizing_fu.volume = V_fu ;                                                     % Save in output structure
sizing_fu.radius = r_tank_fu ;                                                % Save in output structure
sizing_fu.thickness = t_tank_fu ;                                             % Save in output structure
sizing_fu.tank_mass = m_tank_fu;                                              % Save in output structure

sizing_gas.type = pressurant.type ;                                                 % Save in output structure
sizing_gas.mass = m_gas ;                                                     % Save in output structure
sizing_gas.volume = V_gas ;                                                   % Save in output structure
sizing_gas.radius = r_tank_gas ;                                              % Save in output structure
sizing_gas.thickness = t_tank_gas ;                                           % Save in output structure
sizing_gas.tank_mass = m_tank_gas;                                            % Save in output structure

%% Compute the total mass
m_ME = thruster.mass_ME ;                                                      % [ kg ] - Mass of 1 Main Engine
n_ME = thruster.number_ME ;                                                    % [ - ] - Number of Main Engines
m_ST = thruster.mass_ST ;                                                      % [ kg ] - Mass of 1 Small Thruster
n_ST = thruster.number_ST ;                                                    % [ - ] - Number of Small Thrusters
M_PS = m_ox + m_fu + m_gas + m_tank + n_ME * m_ME + n_ST * m_ST ;

end
