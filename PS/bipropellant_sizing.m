function [ M_PS, sizing_ox, sizing_fu, sizing_gas ] = bipropellant_sizing( oxidizer, fuel, pressurant, tank, thruster )

% INPUTS
%
% oxidizer      [ struc ] - Structure with all propellant properties:
%                               oxidizer.rho : density of propellant              [ kg/m^3 ]          
%                               oxidizer.OFratio : O/F ratio                      [ - ]
%                               oxidizer.prop_mass : mass of propellant           [ kg ]
%
% fuel          [ struc ] - Structure with all propellant properties:
%                               fuel.rho : density of propellant                  [ kg/m^3 ]   
%
% pressurant    [ struc ] - Structure with all pressurant properties:
%                               pressurant.feed_system : feeding system type      [ str ]
%                                       - 'Blowdown' 
%                                       - 'PressureRegulated' 
%                               pressurant.feed_system.B : Blowdown ratio (only   [ - ]
%                                                          if blowdown system)
%                               pressurant.gamma : specific heat ratio            [ - ]
%                               pressurant.R : specific gas constant              [ J/kgK ]
%                               pressurant.pressure : initial pressure of         [ Pa ]
%                                                     pressurant     
%
% tank          [ struc ] - Structure with all tank properties
%                               tank.rho : tank material density                  [ kg/m^3 ]
%                               tank.sigma : tank material strenght               [ Pa ]
%                               tank.height : cylindrical intersection height     [ m ]
%                               tank.temperature : tank temperature (if not given [ K ]
%                                                  will be considered as 10 deg)               
%                               tank.type : tank geometry                         [ str ]
%                                       - 'Spherical' 
%                                       - 'Cylindrical' -> cylindrical with
%                                                          spherical caps
%                               tank.number : number of tanks ( defines if        [ - ]
%                                             propellant and pressurant are in 
%                                             the same tank or in saperate tanks)
%                                       - 2: 1 for ox and pressurant, 1 for
%                                            fu and pressurant
%                                       - 3: 1 for ox, 1 for fu and 1 for pressurant
%                                       - 4: 1 for ox, 1 for fu, 1 for oxidizer 
%                                            pressurant and 1 for fuel pressurant
%
% thruster      [ struc ] - Structure with all tank properties
%                               thruster.mass : thruster mass                     [ kg ]
%                               thruster.number : number of thrusters             [ - ]
%                               thruster.chamber_pressure : combustion            [ Pa ]
%                                                           chamber pressure

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

%% Recover data from structures
m_prop = oxidizer.prop_mass ;                                                     % [ kg ] - Propellant mass
m_prop = m_prop * 1.13 ;                                                          % [ kg ] Margin: MAR-CP-010 + 3% of unsable mass
  
OF = oxidizer.OFratio ;                                                           % [ - ] - Oxidizer-to-fuel ratio

m_fu = m_prop / ( 1 + OF ) ;                                                      % [ kg ] - Fuel mass
m_fu = 1.01 * m_fu ;                                                              % [ kg ] - Margin: ullage

m_ox = ( OF * m_prop ) / ( OF + 1 ) ;                                             % [ kg ] - Oxidizer mass
m_ox = 1.01 * m_ox ;                                                              % [ kg ] - Margin: ullage


rho_ox = oxidizer.rho ;                                                           % [ kg/m^3 ] - Oxidizer density
rho_fu = fuel.rho ;                                                               % [ kg/m^3 ] - Fuelt density

rho_prop = m_prop / ( m_fu / rho_fu + m_ox / rho_ox ) ;                           % [ kg/m^3 ] - Mean propellant density

V_ox = m_ox / rho_ox ;                                                            % [ m^3 ] - Oxidizer volume
V_fu = m_fu / rho_fu ;                                                            % [ m^3 ] - Fuel volume

V_prop = m_prop / rho_prop ;                                                      % [ m^3 ] - Propellant volume

% If the initial pressure of the propellant is not given as input, it is
% retrieved from pressure losses and desired combustion chamber pressure 
if isfield( thruster, 'chamber_pressure' )
    P_c = thruster.chamber_pressure ;                                             % [ Pa ] - Combustion chamber pressure
else
    P_c = 20e5 ;        

end

dP_inj = 0.3 * P_c ;                                                              % [ Pa ] - Injection losses
dP_feed = 40e3 ;                                                                  % [ Pa ] - Feeding system losses
Pi_prop = P_c + dP_inj + dP_feed ;                                                % [ Pa ] - Propellant tank initial pressure

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

m_thruster = thruster.mass ;                                                      % [ kg ] - Mass of single thruster
n_thruster = thruster.number ;                                                    % [ - ] - Number of thrustres 

%% Compute gas volume

switch pressurant.feed_system

    case 'Blowdown'
        B = pressurant.feed_system. B ;                                           % [ - ] - Blowdown ratio

        % Oxidizer
        V_gas_ox = V_ox / ( B - 1) ;                                              % [ m^3 ] - Pressurant volume for oxidizer
        V_gas_ox = V_gas_ox * 1.01 ;                                              % [ m^3 ] - Margin: baldder volume

        m_gas_ox = ( Pi_gas * V_gas_ox ) / ( R * T_tank ) ;                       % [ kg ] - Pressurant mass for oxidizer
        m_gas_ox = m_press * 1.2 ;                                                % [ kg ] - Margin: MAR-MAS-090

        % Fuel
        V_gas_fu = V_fu / ( B - 1) ;                                              % [ m^3 ] - Pressurant volume for fuel
        V_gas_fu = V_gas_fu * 1.01 ;                                              % [ m^3 ] - Margin: baldder volume

        m_gas_fu = ( Pi_gas * V_gas_fu ) / ( R * T_tank ) ;                       % [ kg ] - Pressurant mass for fuel
        m_gas_fu = m_press * 1.2 ;                                                % [ kg ] - Margin: MAR-MAS-090

        m_gas = m_gas_ox + m_gas_fu ;

    case 'PressureRegulated'
        % Oxidizer
        m_gas_ox = ( Pi_prop * V_ox ) / ( R * T_tank ) * gamma / ...               % [ kg ] - Pressurant mass
                ( 1 - Pi_prop / Pi_gas ) ;                             
        m_gas_ox = m_gas_ox * 1.2 ;                                                   % [ kg ] - Margin: MAR-MAS-090      

        V_gas_ox = ( m_gas_ox * R * T_tank ) / Pi_gas ;                                 % [ m^3 ] - Pressurant volume (coincides with initial pressurant tank volume)
        V_gas_ox = 1.01 * V_gas_ox ;                                                    % [ m^3 ] - Margin: bladder 

        % Fuel
        m_gas_fu = ( Pi_prop * V_fu ) / ( R * T_tank ) * gamma / ...               % [ kg ] - Pressurant mass
                ( 1 - Pi_prop / Pi_gas ) ;                             
        m_gas_fu = m_gas_fu * 1.2 ;                                                   % [ kg ] - Margin: MAR-MAS-090      

        V_gas_fu = ( m_gas_fu * R * T_tank ) / Pi_gas ;                                 % [ m^3 ] - Pressurant volume (coincides with initial pressurant tank volume)
        V_gas_fu = 1.01 * V_gas_fu ;                                                    % [ m^3 ] - Margin: bladder 

        m_gas = m_gas_ox + m_gas_fu ;
        V_gas = V_gas_ox + V_gas_fu ;
end

%% Compute tank geometry and mass

if tank.number == 2 % Case in which pressurant and propellant are in the same tank (1 for ox, 1 for fu)
    
    switch tank.type

        case 'Spherical'
            % Oxidizer & pressurant
            r_tank_ox = ( ( 3 * ( V_gas_ox + V_ox ) ) / ( 4 * pi ) ) ^ ( 1 / 3 ); % [ m ] - Oxidizer-pressurant tank radius
            t_tank_ox = ( Pi_gas * r_tank_ox ) / ( 2 * sigma_tank ) ;             % [ m ] - Oxidizer-pressurant tank thickness
    
            m_tank_ox = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_ox + ...           % [ kg ] - Oxidizer-pressurant tank mass
                        t_tank_ox ) ^ 3 - r_tank_ox ^ 3 ) ; 

            % Fuel & pressurant
            r_tank_fu = ( ( 3 * ( V_gas_fu + V_fu ) ) / ( 4 * pi ) ) ^ ( 1 / 3 ); % [ m ] - Fuel-pressurant tank radius
            t_tank_fu = ( Pi_gas * r_tank_fu ) / ( 2 * sigma_tank ) ;             % [ m ] - Fuel-pressurant tank thickness
    
            m_tank_fu = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_fu + ...           % [ kg ] - Fuel-pressurant tank mass
                        t_tank_fu ) ^ 3 - r_tank_fu ^ 3 ) ; 
            
            %
            m_tank = m_tank_ox + m_tank_fu ;                                      % [ kg ] - Total tank mass

        case 'Cylindrical'
            h_tank = tank.heigh ;                                                 % [ m ] - Height of cylindrical interection 
            
            % Oxidizer & pressurant
            r_function_ox = @( x ) ( V_gas_ox + V_ox ) - ( 4 * pi / 3 ) * ...     % Function used to compute the radius of the tank knowing the height
                                   x ^ 3 - pi * h_tank * x ^ 2 ; 
            r_guess_ox = 0.3 * h_tank ;
            r_tank_ox = fzero( r_function_ox, r_guess_ox ) ;   

            t_tank_ox = ( Pi_gas * r_tank_ox ) / sigma_tank ;                     % [ m ] - Oxidizer-pressurant tank thickness
           
            m_tank_ox = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_ox +  ...        % [ kg ] - Oxidizer-pressurant tank mass
                        t_tank_ox ) ^ 3 - r_tank_ox ^ 3 ) +  pi * h_tank * ( ...
                        ( r_tank_ox + t_tank_ox) ^ 2 - r_tank_ox ^ 2 ) ) ;

            % Fuel & pressurant
            r_function_fu = @( x ) ( V_gas_fu + V_fu ) - ( 4 * pi / 3 ) *  ...    % Function used to compute the radius of the tank knowing the height
                                x ^ 3 - pi * h_tank * x ^ 2 ; 
            r_guess_fu = 0.3 * h_tank ;
            r_tank_fu = fzero( r_function_fu, r_guess_fu ) ;   

            t_tank_fu = ( Pi_gas * r_tank_fu ) / sigma_tank ;                     % [ m ] - Fuel-pressurant tank thickness
           
            m_tank_fu = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_fu + ...         % [ kg ] - Fuel-pressurant tank mass
                        t_tank_fu ) ^ 3 - r_tank_fu ^ 3 ) +  pi * h_tank * ( ...
                        ( r_tank_fu + t_tank_fu) ^ 2 - r_tank_fu ^ 2 ) ) ;

            %
            m_tank = m_tank_ox + m_tank_fu ;                                      % [ kg ] - Total tank mass
            
    end


    sizing_ox.single.type = tank.type ;                                           % Save in output structure
    sizing_ox.single.mass_press = m_gas_ox ;                                      % Save in output structure
    sizing_ox.single.volume_press = V_gas_ox ;                                    % Save in output structure
    sizing_ox.single.mass_ox = m_ox ;                                             % Save in output structure
    sizing_ox.single.volume_ox = V_ox ;                                           % Save in output structure
    sizing_ox.single.radius = r_tank_ox ;                                         % Save in output structure
    sizing_ox.single.thickness = t_tank_ox ;                                      % Save in output structure
    sizing_ox.single.tank_mass = m_tank_ox ;                                      % Save in output structure
    sizing_ox.single.tank_volume = V_gas_ox + V_ox ;                              % Save in output structure

    sizing_fu.single.type = tank.type ;                                           % Save in output structure
    sizing_fu.single.mass_press = m_gas_fu ;                                      % Save in output structure
    sizing_fu.single.volume_press = V_gas_fu ;                                    % Save in output structure
    sizing_fu.single.mass_fu = m_fu ;                                             % Save in output structure
    sizing_fu.single.volume_fu = V_fu ;                                           % Save in output structure
    sizing_fu.single.radius = r_tank_fu ;                                         % Save in output structure
    sizing_fu.single.thickness = t_tank_fu ;                                      % Save in output structure
    sizing_fu.single.tank_mass = m_tank_fu ;                                      % Save in output structure
    sizing_fu.single.tank_volume = V_gas_fu + V_fu ;                              % Save in output structure

elseif tank.number == 3  % Case in which pressurant and propellant are in separate tanks ( 1 common pressurant, 1 ox and 1 fu ) 

    switch tank.type

        case 'Spherical'
            % Pressurant tank
            r_tank_gas = ( ( 3 * ( V_gas_fu + V_gas_ox ) ) / ( 4 * pi ) ) ^ ...   % [ m ] - Pressurant tank radius
                         ( 1 / 3 );             
            t_tank_gas = ( Pi_gas * r_tank_gas ) / ( 2 * sigma_tank ) ;           % [ m ] - Pressurant anks thickness
            m_tank_gas = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_gas + ...         % [ kg ] - Pressurant tank mass 
                         t_tank_gas ) ^ 3 - r_tank_gas ^ 3 ) ; 

            % Oxidizer tank
            r_tank_ox = ( ( 3 * V_ox ) / ( 4 * pi ) ) ^ ( 1 / 3 );                % [ m ] - Propellant tank radius
            t_tank_ox = ( Pi_prop * r_tank_ox ) / ( 2 * sigma_tank ) ;            % [ m ] - Propellant tank thickness
            m_tank_ox = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_ox + ...           % [ kg ] - Propellant tank mass
                          t_tank_ox) ^ 3 - r_tank_ox ^ 3 ) ;

            % Fuel tank
            r_tank_fu = ( ( 3 * V_fu ) / ( 4 * pi ) ) ^ ( 1 / 3 );                % [ m ] - Propellant tank radius
            t_tank_fu = ( Pi_prop * r_tank_fu ) / ( 2 * sigma_tank ) ;            % [ m ] - Propellant tank thickness
            m_tank_fu = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_fu + ...           % [ kg ] - Propellant tank mass
                          t_tank_fu) ^ 3 - r_tank_fu ^ 3 ) ;

            m_tank = m_tank_gas + m_tank_ox + m_tank_fu ;                         % [ kg ] - Total tank mass

        case 'Cylindrical'
    
            h_tank = tank.heigh ;                                                 % [ m ] - Height of cylindrical interection - NOTE: we assumed both propellant and pressurant tanks to have same height
            
            % Pressurant tank
            r_function_gas = @( x ) ( V_gas_ox + V_gas_fu ) - ( 4 * pi / 3 )  ... % Function used to compute the radius of the tank knowing the height
                                    * x ^ 3 - pi * h_tank * x ^ 2 ; 
            r_guess_gas = 0.3 * h_tank ;
            r_tank_gas = fzero( r_function_gas, r_guess_gas ) ;   
            t_tank_gas = ( Pi_gas * r_tank_gas ) / sigma_tank ;                   % [ m ] - Tank thickness
           
            m_tank_gas = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_gas + ...       % [ kg ] - Tank mass
                         t_tank_gas ) ^ 3 - r_tank_gas ^ 3 ) +  pi * h_tank *...
                         ( ( r_tank_gas + t_tank_gas) ^ 2 - r_tank_gas ^ 2 ) ) ;

            % Oxidizer tank
            r_function_ox = @( x ) V_ox - ( 4 * pi / 3 ) * x ^ 3 - pi * ...       % Function used to compute the radius of the tank knowing the height
                                     h_tank * x ^ 2 ; 
            r_guess_ox = 0.3 * h_tank ;
            r_tank_ox = fzero( r_function_ox, r_guess_ox ) ;   

            t_tank_ox = ( Pi_prop * r_tank_ox ) / sigma_tank ;                    % [ m ] - Tank thickness
           
            m_tank_ox = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_ox + ...         % [ kg ] - Tank mass
                         t_tank_ox ) ^ 3 - r_tank_ox ^ 3 ) +  pi * h_tank * ...
                         ( ( r_tank_ox + t_tank_ox ) ^ 2 - r_tank_ox ^ 2 ) ) ;
 
            % Fuel tank
            r_function_fu = @( x ) V_fu - ( 4 * pi / 3 ) * x ^ 3 - pi * ...       % Function used to compute the radius of the tank knowing the height
                                     h_tank * x ^ 2 ; 
            r_guess_fu = 0.3 * h_tank ;
            r_tank_fu = fzero( r_function_fu, r_guess_fu ) ;   

            t_tank_fu = ( Pi_prop * r_tank_fu ) / sigma_tank ;                    % [ m ] - Tank thickness
           
            m_tank_fu = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_fu + ...         % [ kg ] - Tank mass
                         t_tank_fu ) ^ 3 - r_tank_fu ^ 3 ) +  pi * h_tank * ...
                         ( ( r_tank_fu + t_tank_fu ) ^ 2 - r_tank_fu ^ 2 ) ) ;

            %
            m_tank = m_tank_gas + m_tank_ox + m_tank_fu ;                         % [ kg ] - Total tank mass

    end

    sizing_ox.type = tank.type ;                                                  % Save in output structure
    sizing_ox.mass = m_ox ;                                                       % Save in output structure
    sizing_ox.volume = V_ox ;                                                     % Save in output structure
    sizing_ox.radius = r_tank_ox ;                                                % Save in output structure
    sizing_ox.thickness = t_tank_ox ;                                             % Save in output structure
    sizing_ox.tank_mass = m_tank_ox ;                                             % Save in output structure

    sizing_fu.type = tank.type ;                                                  % Save in output structure
    sizing_fu.mass = m_fu ;                                                       % Save in output structure
    sizing_fu.volume = V_fu ;                                                     % Save in output structure
    sizing_fu.radius = r_tank_fu ;                                                % Save in output structure
    sizing_fu.thickness = t_tank_fu ;                                             % Save in output structure
    sizing_fu.tank_mass = m_tank_fu;                                              % Save in output structure

    sizing_gas.type = tank.type ;                                                 % Save in output structure
    sizing_gas.mass = m_gas ;                                                     % Save in output structure
    sizing_gas.volume = V_gas ;                                                   % Save in output structure
    sizing_gas.radius = r_tank_gas ;                                              % Save in output structure
    sizing_gas.thickness = t_tank_gas ;                                           % Save in output structure
    sizing_gas.tank_mass = m_tank_gas;                                            % Save in output structure

    sizing.propellant.total_mass = m_tank ;                                       % Save in output structure

elseif tank.number == 4  % Case in which pressurant and propellant are in separate tanks ( 1 pressurant for ox, 1 pressurant for ru, 1 ox and 1 fu ) 

    switch tank.type

        case 'Spherical'
            % Oxider pressurant tank
            r_tank_gas_ox = ( ( 3 * V_gas_ox ) / ( 4 * pi ) ) ^ ...               % [ m ] - Pressurant tank radius
                         ( 1 / 3 );             
            t_tank_gas_ox = ( Pi_gas * r_tank_gas_ox ) / ( 2 * sigma_tank ) ;     % [ m ] - Pressurant anks thickness
            m_tank_gas_ox = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_gas_ox + ...   % [ kg ] - Pressurant tank mass 
                         t_tank_gas_ox ) ^ 3 - r_tank_gas_ox ^ 3 ) ; 

            % Oxidizer tank
            r_tank_ox = ( ( 3 * V_ox ) / ( 4 * pi ) ) ^ ( 1 / 3 );                % [ m ] - Propellant tank radius
            t_tank_ox = ( Pi_prop * r_tank_ox ) / ( 2 * sigma_tank ) ;            % [ m ] - Propellant tank thickness
            m_tank_ox = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_ox + ...           % [ kg ] - Propellant tank mass
                          t_tank_ox) ^ 3 - r_tank_ox ^ 3 ) ;

            % Fuel pressurant tank
            r_tank_gas_fu = ( ( 3 * V_gas_fu ) / ( 4 * pi ) ) ^ ...               % [ m ] - Pressurant tank radius
                         ( 1 / 3 );             
            t_tank_gas_fu = ( Pi_gas * r_tank_gas_fu ) / ( 2 * sigma_tank ) ;     % [ m ] - Pressurant anks thickness
            m_tank_gas_fu = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_gas_fu + ...   % [ kg ] - Pressurant tank mass 
                         t_tank_gas_fu ) ^ 3 - r_tank_gas_fu ^ 3 ) ; 

            % Fuel tank
            r_tank_fu = ( ( 3 * V_fu ) / ( 4 * pi ) ) ^ ( 1 / 3 );                % [ m ] - Propellant tank radius
            t_tank_fu = ( Pi_prop * r_tank_fu ) / ( 2 * sigma_tank ) ;            % [ m ] - Propellant tank thickness
            m_tank_fu = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_fu + ...           % [ kg ] - Propellant tank mass
                          t_tank_fu) ^ 3 - r_tank_fu ^ 3 ) ;

            m_tank = m_tank_gas_ox + m_tank_gas_fu + m_tank_ox + m_tank_fu ;      % [ kg ] - Total tank mass

        case 'Cylindrical'
    
            h_tank = tank.heigh ;                                                 % [ m ] - Height of cylindrical interection - NOTE: we assumed both propellant and pressurant tanks to have same height

            % Oxidizer pressurant tank
            r_function_gas_ox = @( x ) V_gas_ox - ( 4 * pi / 3 ) * x ^ 3 - ...    % Function used to compute the radius of the tank knowing the height
                                    pi * h_tank * x ^ 2 ; 
            r_guess_gas_ox = 0.3 * h_tank ;
            r_tank_gas_ox = fzero( r_function_gas_ox, r_guess_gas_ox ) ;   
            t_tank_gas_ox = ( Pi_gas * r_tank_gas_ox ) / sigma_tank ;             % [ m ] - Oxidizer pressurant tank thickness
           
            m_tank_gas_ox = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_gas_ox + ... % [ kg ] - Oxidizer pressurant tank mass
                         t_tank_gas_ox ) ^ 3 - r_tank_gas_ox ^ 3 ) +  pi * h_tank *...
                         ( ( r_tank_gas_ox + t_tank_gas_ox ) ^ 2 - r_tank_gas_ox ^ 2 ) ) ;

            % Oxidizer tank
            r_function_ox = @( x ) V_ox - ( 4 * pi / 3 ) * x ^ 3 - pi * ...       % Function used to compute the radius of the tank knowing the height
                                     h_tank * x ^ 2 ; 
            r_guess_ox = 0.3 * h_tank ;
            r_tank_ox = fzero( r_function_ox, r_guess_ox ) ;   

            t_tank_ox = ( Pi_prop * r_tank_ox ) / sigma_tank ;                    % [ m ] - Oxidizer tank thickness
           
            m_tank_ox = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_ox + ...         % [ kg ] - Oxidizer tank mass
                         t_tank_ox ) ^ 3 - r_tank_ox ^ 3 ) +  pi * h_tank * ...
                         ( ( r_tank_ox + t_tank_ox ) ^ 2 - r_tank_ox ^ 2 ) ) ;
 
            % Fuel pressurant tank
            r_function_gas_fu = @( x ) V_gas_fu - ( 4 * pi / 3 ) * x ^ 3 - ...    % Function used to compute the radius of the tank knowing the height
                                    pi * h_tank * x ^ 2 ; 
            r_guess_gas_fu = 0.3 * h_tank ;
            r_tank_gas_fu = fzero( r_function_gas_fu, r_guess_gas_fu ) ;   
            t_tank_gas_fu = ( Pi_gas * r_tank_gas_fu ) / sigma_tank ;             % [ m ] - Fuel pressurant tank thickness
           
            m_tank_gas_fu = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_gas_fu + ... % [ kg ] - Fuel pressurant tank mass
                         t_tank_gas_fu ) ^ 3 - r_tank_gas_fu ^ 3 ) +  pi * h_tank *...
                         ( ( r_tank_gas_fu + t_tank_gas_fu ) ^ 2 - r_tank_gas_fu ^ 2 ) ) ;

            % Fuel tank
            r_function_fu = @( x ) V_fu - ( 4 * pi / 3 ) * x ^ 3 - pi * ...       % Function used to compute the radius of the tank knowing the height
                                     h_tank * x ^ 2 ; 
            r_guess_fu = 0.3 * h_tank ;
            r_tank_fu = fzero( r_function_fu, r_guess_fu ) ;   

            t_tank_fu = ( Pi_prop * r_tank_fu ) / sigma_tank ;                    % [ m ] - Fuel tank thickness
           
            m_tank_fu = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_fu + ...         % [ kg ] - Fuel tank mass
                         t_tank_fu ) ^ 3 - r_tank_fu ^ 3 ) +  pi * h_tank * ...
                         ( ( r_tank_fu + t_tank_fu ) ^ 2 - r_tank_fu ^ 2 ) ) ;

            %
            m_tank = m_tank_gas_ox + m_tank_gas_fu + m_tank_ox + m_tank_fu ;      % [ kg ] - Total tank mass

    end 

    sizing_ox.type = tank.type ;                                                  % Save in output structure
    sizing_ox.mass_press = m_gas_ox ;                                             % Save in output structure
    sizing_ox.volume_press = V_gas_ox ;                                           % Save in output structure
    sizing_ox.radius_press = r_tank_gas_ox ;                                      % Save in output structure
    sizing_ox.thickness_press = t_tank_gas_ox ;                                   % Save in output structure
    sizing_ox.tank_mass_press = m_tank_gas_ox ;                                   % Save in output structure

    sizing_ox.mass = m_ox ;                                                    % Save in output structure
    sizing_ox.volume = V_ox ;                                                  % Save in output structure                                  % Save in output structure
    sizing_ox.radius = r_tank_gas_ox ;                                         % Save in output structure
    sizing_ox.thickness = t_tank_gas_ox ;                                      % Save in output structure
    sizing_ox.tank_mass = m_tank_gas_ox ;                                      % Save in output structure

    sizing_fu.type = tank.type ;                                                  % Save in output structure
    sizing_fu.mass_press = m_gas_fu ;                                             % Save in output structure
    sizing_fu.volume_press = V_gas_fu ;                                           % Save in output structure
    sizing_fu.radius_press = r_tank_gas_fu ;                                      % Save in output structure
    sizing_fu.thickness_press = t_tank_gas_fu ;                                   % Save in output structure
    sizing_fu.tank_mass_press = m_tank_gas_fu ;                                   % Save in output structure

    sizing_fu.mass = m_fu ;                                                    % Save in output structure
    sizing_fu.volume = V_fu ;                                                  % Save in output structure                                  % Save in output structure
    sizing_fu.radius = r_tank_gas_ox ;                                         % Save in output structure
    sizing_fu.thickness = t_tank_gas_fu ;                                      % Save in output structure
    sizing_fu.tank_mass = m_tank_gas_fu ;                                      % Save in output structure

end
%% Compute the total mass
M_PS = m_prop + m_gas + m_tank + n_thruster * m_thruster ;

end