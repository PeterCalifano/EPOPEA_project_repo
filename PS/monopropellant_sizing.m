function [ sizing, M_PS ] = monopropellant_sizing( propellant, pressurant, tank, thruster )

% INPUTS
%
% propellant    [ struc ] - Structure with all propellant properties:
%                               propellant.mass : mass of propellant              [ kg ]
%                               propellant.rho : density of propellant            [ kg/m^3 ]  
%                               propellant.type : geometry of  propellant
%                                                 tank ('Spherical' or 'Cylindrical')
%
% pressurant    [ struc ] - Structure with all pressurant properties:
%                               pressurant.feed_system : feeding system type      [ str ]
%                                       - 'Blowdown' 
%                                       - 'PressureRegulated' 
%                               pressurant.B : Blowdown ratio (only   [ - ]
%                                                          if blowdown system)
%                               pressurant.gamma : specific heat ratio            [ - ]
%                               pressurant.R : specific gas constant              [ J/kgK ]
%                               pressurant.pressure : initial pressure of         [ Pa ]
%                                                     pressurant     
%                               pressurant.type : geometry of pressurant
%                                                 tank ('Spherical' or
%                                                 'Cylindrical') ( only if separate tanks)
%
% tank          [ struc ] - Structure with all tank properties
%                               tank.rho : tank material density                  [ kg/m^3 ]
%                               tank.sigma : tank material strenght               [ MPa ]
%                               tank.temperature : tank temperature (if not given [ K ]
%                                                  will be considered as 10 deg)       
%                               tank.ratio : cylindrical intersection radius-to-height ratio     [ m ]
%                               tank.type : tank geometry                         [ str ]
%                                       - 'Spherical' 
%                                       - 'Cylindrical' -> cylindrical with
%                                                          spherical caps
%                               tank.number : number of tanks ( defines if        [ - ]
%                                             propellant and pressurant are in 
%                                             the same tank or in saperate tanks) 
%
% thruster      [ struc ] - Structure with all tank properties
%                               thruster.mass : thruster mass                     [ kg ]
%                               thruster.number : number of thrusters             [ - ]
%                               thruster.chamber_pressure : combustion            [ Pa ]
%                                                           chamber pressure
%                               thurster.chamber_pressure_min : maximum           [ Pa ]
%                                                               chamber pressure

% OUTPUTS
%
% sizing        [ struc ] - Structure with tank sizing:
%                               sizing.single : tank properties (when single tank [ struc ] 
%                                               is used)         
%                                   - .type                                       [ str ]
%                                   - .mass_press                                 [ kg ]
%                                   - .volume_press                               [ m^3 ]
%                                   - .mass_prop                                  [ kg ]
%                                   - .volume_prop                                [ m^3 ]
%                                   - .radius                                     [ m ]
%                                   - .thickness                                  [ m ]
%                                   - .tank_mass                                  [ kg ]
%                                   - .tank_volume                                [ m^3 ]

%                               sizing.pressurant : pressurant tank properties    [ struc ]
%                                                  (when different tanks are used)
%                                   - .type                                       [ str ]
%                                   - .mass                                       [ kg ]
%                                   - .volume                                     [ m^3 ]
%                                   - .radius                                     [ m ]
%                                   - .thickness                                  [ m ]
%                                   - .tank_mass                                  [ kg ]

%                               sizing.propellant : propellant tank properties    [ struc ]
%                                                   (when different tanks are used)
%                                   - .type                                       [ str ]
%                                   - .mass                                       [ kg ]
%                                   - .volume                                     [ m^3 ]
%                                   - .radius                                     [ m ]
%                                   - .thickness                                  [ m ]
%                                   - .tank_mass                                  [ kg ]
%                                   - .total_mass                                 [ kg ]

%                               NOTE : The total tank mass (when different tanks are used) is saved as
%                                      sizing.propellant.total_mass
%
% M_PS          [ 1 ]     - Total mass of the propulsion subsystem                [ kg ]
% ---

%% Recover data from structures
m_prop = propellant.mass ;                                                        % [ kg ] - Propellant mass
m_prop = m_prop * 1.12 ;                                                          % [ kg ] Margin: MAR-CP-010 + 3% of unsable mass
rho_prop = propellant.rho ;                                                       % [ kg/m^3 ] - Propellant density
V_prop = m_prop / rho_prop ;                                                      % [ m^3 ] - Propellant volume

% The initial pressure of the propellant is retrieved from pressure losses and
% desired combustion chamber pressure 
if isfield( thruster, 'chamber_pressure' )
    P_c = thruster.chamber_pressure ;                                             % [ Pa ] - Combustion chamber pressure
else
    P_c = 20e5 ;        

end

dP_inj = 0.3 * P_c ;                                                              % [ Pa ] - Injection losses
dP_feed = 50e3 ;                                                                  % [ Pa ] - Feeding system losses
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
        B = pressurant.B ;                                                        % [ - ] - Blowdown ratio
        
        P_c_min = thruster.chamber_pressure_min ;
        Pf_gas = P_c_min + dP_feed ;
        B_max = Pi_gas / Pf_gas ;

        if B > B_max
            error( 'Blowdown ratio exceeds the maximum value: %d \n', B_max ) 
        end

        V_gas = V_prop / ( B - 1) ;                                            % [ m^3 ] - Pressurant volume
        V_gas = V_gas * 1.01 ;                                                    % [ m^3 ] - Margin: baldder volume

        m_gas = ( Pi_gas * V_gas ) / ( R * T_tank ) ;                             % [ kg ] - Pressurant mass
        m_gas = m_gas * 1.2 ;                                                   % [ kg ] - Margin: MAR-MAS-090

    case 'PressureRegulated'
        Pf_gas = 1 / 10 * Pi_gas ;
        m_gas = ( Pi_prop * V_prop ) / ( R * T_tank ) * gamma / ...               % [ kg ] - Pressurant mass
                ( 1 - Pf_gas / Pi_gas ) ;                             
        m_gas = m_gas * 1.2 ;                                                   % [ kg ] - Margin: MAR-MAS-090      

        V_gas = ( m_gas * R * T_tank ) / Pi_gas ;                                 % [ m^3 ] - Pressurant volume (coincides with initial pressurant tank volume)
        V_gas = 1.01 * V_gas ;                                                    % [ m^3 ] - Margin: bladder 

end

%% Compute tank geometry and mass

if tank.number == 1 % Case in which pressurant and propellant are in the same tank 
    
    switch propellant.type

        case 'Spherical'
            r_tank = ( ( 3 * ( V_gas + V_prop ) ) / ( 4 * pi ) ) ^ ( 1 / 3 );     % [ m ] - Tank radius
            t_tank = ( Pi_gas * r_tank ) / ( 2 * sigma_tank ) ;                   % [ m ] - Tank thickness
    
            m_tank = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank + t_tank ) ^ 3 - ...  % [ kg ] - Tank mass
                     r_tank ^ 3 ) ; 

        case 'Cylindrical'
            ratio = tank.ratio ;                                                 % [ m ] - Radius-to-height ratio of cylindrical interection - NOTE: we assumed both propellant and pressurant tanks to have same height
            % height = radius / ratio

            r_function = @( x ) ( V_gas + V_prop ) - ( 4 * pi / 3 ) * x ^ 3 - ... % Function used to compute the radius of the tank knowing the height
                                pi * x / ratio * x ^ 2 ; 
            r_guess = 0.3 * ratio ;
            r_tank = fzero( r_function, r_guess ) ;   
            h_tank = r_tank / ratio ;
            sizing.heigth = h_tank ;                                    % Save in output structure
            
            t_tank = ( Pi_gas * r_tank ) / sigma_tank ;                           % [ m ] - Tank thickness
           
            m_tank = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank + t_tank ) ^ 3 - ...% [ kg ] - Tank mass
                     r_tank ^ 3 ) +  pi * h_tank * ( ( r_tank + t_tank) ^ 2 - r_tank ^ 2 ) ) ;

    end

    sizing.type = propellant.type ;                                              % Save in output structure
    sizing.mass_press = m_gas ;                                            % Save in output structure
    sizing.volume_press = V_gas ;                                          % Save in output structure
    sizing.mass_prop = m_prop ;                                            % Save in output structure
    sizing.volume_prop = V_prop ;                                          % Save in output structure
    sizing.radius = r_tank ;                                               % Save in output structure
    sizing.thickness = t_tank ;                                            % Save in output structure
    sizing.tank_mass = m_tank ;                                            % Save in output structure
    sizing.tank_volume = V_gas + V_prop ;                                  % Save in output structure


%     switch tank.type
% 
%         case 'Spherical'
%             r_tank = ( ( 3 * ( V_gas + V_prop ) ) / ( 4 * pi ) ) ^ ( 1 / 3 );     % [ m ] - Tank radius
%             t_tank = ( Pi_gas * r_tank ) / ( 2 * sigma_tank ) ;                   % [ m ] - Tank thickness
%     
%             m_tank = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank + t_tank ) ^ 3 - ...  % [ kg ] - Tank mass
%                      r_tank ^ 3 ) ; 

%         case 'Cylindrical'
%             ratio = tank.ratio ;                                                 % [ m ] - Radius-to-height ratio of cylindrical interection - NOTE: we assumed both propellant and pressurant tanks to have same height
% 
%             % height = radius / ratio
% 
%             r_function = @( x ) ( V_gas + V_prop ) - ( 4 * pi / 3 ) * x ^ 3 - ... % Function used to compute the radius of the tank knowing the height
%                                 pi * x / ratio * x ^ 2 ; 
%             r_guess = 0.3 * ratio ;
%             r_tank = fzero( r_function, r_guess ) ;   
%             h_tank = r_tank / ratio ;
%             t_tank = ( Pi_gas * r_tank ) / sigma_tank ;                           % [ m ] - Tank thickness
%            
%             m_tank = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank + t_tank ) ^ 3 - ...% [ kg ] - Tank mass
%                      r_tank ^ 3 ) +  pi * h_tank * ( ( r_tank + t_tank) ^ 2 - r_tank ^ 2 ) ) ;
%             sizing.single.heigth = r_tank / ratio ;                                    % Save in output structure
%             
% end

else % Case in which pressurant and propellant are in separate tanks  

    switch propellant.type

        case 'Spherical'
            r_tank_prop = ( ( 3 * V_prop ) / ( 4 * pi ) ) ^ ( 1 / 3 );            % [ m ] - Propellant tank radius
            t_tank_prop = ( Pi_prop * r_tank_prop ) / ( 2 * sigma_tank ) ;        % [ m ] - Propellant tank thickness
            m_tank_prop = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_prop + ...       % [ kg ] - Propellant tank mass
                          t_tank_prop ) ^ 3 - r_tank_prop ^ 3 ) ;

        case 'Cylindrical'
            ratio = tank.ratio ;                                                  % [ m ] - Radius-to-height ratio of cylindrical interection - NOTE: we assumed both propellant and pressurant tanks to have same height
            % height = radius / ratio

            r_function_prop = @( x ) V_prop - ( 4 * pi / 3 ) * x ^ 3 - pi * ...   % Function used to compute the radius of the tank knowing the height
                                     x / ratio * x ^ 2 ; 
            r_guess_prop = 0.3 * ratio ;
            r_tank_prop = fzero( r_function_prop, r_guess_prop) ;   
            h_tank = r_tank_prop / ratio ;
            sizing.pressurant.heigth = h_tank ;                                    % Save in output structure
    
            t_tank_prop = ( Pi_prop * r_tank_prop ) / sigma_tank ;                % [ m ] - Tank thickness
           
            m_tank_prop = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_prop + ...     % [ kg ] - Tank mass
                         t_tank_prop ) ^ 3 - r_tank_prop ^ 3 ) +  pi * h_tank * ...
                         ( ( r_tank_prop + t_tank_prop) ^ 2 - r_tank_prop ^ 2 ) ) ;

    end 

    switch pressurant.type

        case 'Spherical'
            r_tank_gas = ( ( 3 * V_gas ) / ( 4 * pi ) ) ^ ( 1 / 3 );              % [ m ] - Pressurant tank radius
            t_tank_gas = ( Pi_gas * r_tank_gas ) / ( 2 * sigma_tank ) ;           % [ m ] - Pressurant anks thickness
            m_tank_gas = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_gas + ...         % [ kg ] - Pressurant tank mass 
                         t_tank_gas ) ^ 3 - r_tank_gas ^ 3 ) ; 

        case 'Cylindrical'

            ratio = tank.ratio ;                                                  % [ m ] - Radius-to-height ratio of cylindrical interection - NOTE: we assumed both propellant and pressurant tanks to have same height
            % height = radius / ratio

            r_function_gas = @( x ) V_gas - ( 4 * pi / 3 ) * x ^ 3 - ...          % Function used to compute the radius of the tank knowing the height
                                    pi * x / ratio * x ^ 2 ; 
            r_guess_gas = 0.3 * ratio ;
            r_tank_gas = fzero( r_function_gas, r_guess_gas ) ;   
            h_tank = r_tank_gas / ratio ;
            sizing.pressurant.heigth = h_tank ;                                   % Save in output structure
            t_tank_gas = ( Pi_gas * r_tank_gas ) / sigma_tank ;                   % [ m ] - Tank thickness
           
            m_tank_gas = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_gas + ...       % [ kg ] - Tank mass
                         t_tank_gas ) ^ 3 - r_tank_gas ^ 3 ) +  pi * h_tank *...
                         ( ( r_tank_gas + t_tank_gas) ^ 2 - r_tank_gas ^ 2 ) ) ;


    end

    m_tank = m_tank_gas + m_tank_prop ;                                   % [ kg ] - Total tank mass

    sizing.pressurant.type = pressurant.type ;                                          % Save in output structure
    sizing.pressurant.mass = m_gas ;                                              % Save in output structure
    sizing.pressurant.volume = V_gas ;                                            % Save in output structure
    sizing.pressurant.radius = r_tank_gas ;                                       % Save in output structure
    sizing.pressurant.thickness = t_tank_gas ;                                    % Save in output structure
    sizing.pressurant.tank_mass = m_tank_gas ;                                    % Save in output structure

    sizing.propellant.type = propellant.type ;                                          % Save in output structure
    sizing.propellant.mass = m_prop ;                                             % Save in output structure
    sizing.propellant.volume = V_prop ;                                           % Save in output structure
    sizing.propellant.radius = r_tank_prop ;                                      % Save in output structure
    sizing.propellant.thickness = t_tank_prop ;                                   % Save in output structure
    sizing.propellant.tank_mass = m_tank_prop;                                    % Save in output structure


%     switch tank.type

%         case 'Spherical'
            % Pressurant tank
%             r_tank_gas = ( ( 3 * V_gas ) / ( 4 * pi ) ) ^ ( 1 / 3 );              % [ m ] - Pressurant tank radius
%             t_tank_gas = ( Pi_gas * r_tank_gas ) / ( 2 * sigma_tank ) ;           % [ m ] - Pressurant anks thickness
%             m_tank_gas = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_gas + ...         % [ kg ] - Pressurant tank mass 
%                          t_tank_gas ) ^ 3 - r_tank_gas ^ 3 ) ; 

            % Propellant tank
%             r_tank_prop = ( ( 3 * V_prop ) / ( 4 * pi ) ) ^ ( 1 / 3 );            % [ m ] - Propellant tank radius
%             t_tank_prop = ( Pi_prop * r_tank_prop ) / ( 2 * sigma_tank ) ;        % [ m ] - Propellant tank thickness
%             m_tank_prop = rho_tank * ( 4 * pi / 3 ) * ( ( r_tank_prop + ...       % [ kg ] - Propellant tank mass
%                           t_tank_prop ) ^ 3 - r_tank_prop ^ 3 ) ;

%             m_tank = m_tank_gas + m_tank_prop ;                                   % [ kg ] - Total tank mass

%         case 'Cylindrical'
%     
%             ratio = tank.ratio ;                                                 % [ m ] - Radius-to-height ratio of cylindrical interection - NOTE: we assumed both propellant and pressurant tanks to have same height
% 
%             % height = radius / ratio
% 
%             % Pressurant tank
%             r_function_gas = @( x ) V_gas - ( 4 * pi / 3 ) * x ^ 3 - ...          % Function used to compute the radius of the tank knowing the height
%                                     pi * x / ratio * x ^ 2 ; 
%             r_guess_gas = 0.3 * ratio ;
%             r_tank_gas = fzero( r_function_gas, r_guess_gas ) ;   
%             h_tank = r_tank_gas / ratio ;
%             t_tank_gas = ( Pi_gas * r_tank_gas ) / sigma_tank ;                   % [ m ] - Tank thickness
%            
%             m_tank_gas = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_gas + ...       % [ kg ] - Tank mass
%                          t_tank_gas ) ^ 3 - r_tank_gas ^ 3 ) +  pi * h_tank *...
%                          ( ( r_tank_gas + t_tank_gas) ^ 2 - r_tank_gas ^ 2 ) ) ;
% 
%             sizing.pressurant.heigth = r_tank_gas / ratio ;                                    % Save in output structure

            % Propellant tank
%             r_function_prop = @( x ) V_prop - ( 4 * pi / 3 ) * x ^ 3 - pi * ...   % Function used to compute the radius of the tank knowing the height
%                                      h_tank * x ^ 2 ; 
%             r_guess_prop = 0.3 * h_tank ;
%             r_tank_prop = fzero( r_function_prop, r_guess_prop) ;   
%             h_tank = r_tank_gas / ratio ;
%             t_tank_prop = ( Pi_prop * r_tank_prop ) / sigma_tank ;                % [ m ] - Tank thickness
%            
%             m_tank_prop = rho_tank * ( ( 4 * pi / 3 ) * ( ( r_tank_prop + ...     % [ kg ] - Tank mass
%                          t_tank_prop ) ^ 3 - r_tank_prop ^ 3 ) +  pi * h_tank * ...
%                          ( ( r_tank_prop + t_tank_prop) ^ 2 - r_tank_prop ^ 2 ) ) ;
% 
%             m_tank = m_tank_gas + m_tank_prop ;                                   % [ kg ] - Total tank mass
% 
%             sizing.pressurant.heigth = r_tank_prop / ratio ;                                    % Save in output structure
%     end

%     sizing.propellant.total_mass = m_tank ;                                       % Save in output structure

end 

%% Compute the total mass
M_PS = m_prop + m_gas + m_tank + n_thruster * m_thruster ;

% 5% margin shall be added to include pipes
end