function [propellant_tank,pressurant_tank,m_gas] = biprop_sizing(mp,propellant,prop_tank,pressurant,gas_tank)

% Add margins according to CDF standards 
mp = mp*1.13;                                 % Margin MAR-CP-010 + 3% "unusable" propellant / ullage
% Propellant characteristics
rho_fu = propellant.rho_fu;
rho_ox = propellant.rho_ox;
of = propellant.ratio;                        % Take average values from slides (for now)
% Propellant tank properties
Pi_p = 22*1e+5;                               % average from datasheet: 
                                              % PROPELLANT TANK FOR SPINSTABILIZED SPACECRAFT T 11/0
                                              % Check that is > than Pc of thrusters!
sigma_p = prop_tank.sigma;
rho_tank_p = prop_tank.rho;
% Pressurant characteristics
R_gas = pressurant.R; gamma_gas = pressurant.gamma;
% Pressurant tanks initial properties
Pi_gas = 250*1e+5;                            % Change this value! Its an average from datasheets            
% Pressurant tank characteristics
sigma_gas = gas_tank.sigma;
rho_tank_gas = gas_tank.rho;

% Rt_real = gas_tank.radius;
% Vt_real = gas_tank.volume;

% Assumptions
T_tank = 293;                   % Ambient temperature [K]
switch pressurant.mode
    case 'Blowdown'
        % PROPELLANT
        m_fu = (1/(of+1))*mp;                               % fuel mass
        m_ox = (of/(of+1))*mp;                              % oxidizer mass
        rho_mean = mp/(m_fu/rho_fu + m_ox/rho_ox);          % propellant mean density
        V_p = mp/rho_mean;
%         if V_p < 0.3                                     % Average volume limit assumed from Ariane datasheet
            propellant_tank.shape = 'spherical';
            % Spherical
            Rt_prop = (3*V_p/(4*pi))^(1/3);
            t_prop = Pi_p*Rt_prop/(2*sigma_p);
            mt_prop = rho_tank_p*(4*pi/3)*((Rt_prop+t_prop)^3 - Rt_prop^3);
%         else 
%             propellant_tank.shape = 'cylindrical';
%             % Cylindrical with caps
%h_prop = prop_tank.height;        % Different values of the height of cylinder tank, only for this preliminary phase
%             f = @(x) V_p - (4*pi/3)*x^3 - pi*h_prop*x^2;
%             Rt_prop = fzero(@(x) f, 0.3*h_prop);
%             t_prop = Pi_p*Rt_prop/(sigma_p);
%             mt_prop = rho_tank_p*((4*pi/3)*((Rt_prop+t_prop)^3 - Rt_prop^3) +...
%                 pi*h_prop*(Rt_prop+t_prop)^2 - Rt_prop^2);
%         end

        % PRESSURANT
        B = pressurant.B;                                     % Chosen blowdown ratio
        Vt_gas = V_p/(B-1);                                 % initial volume of gas needed
        m_gas = (Pi_gas*Vt_gas)/(R_gas*T_tank);             % Gas initial and tank volume
        m_gas = m_gas*1.2;                                  % Margin MAR-MAS-090
        
%         if Vt_gas < 0.3                                     % Average volume limit assumed from Ariane datasheet
            pressurant_tank.shape = 'spherical';
            % Spherical case
            Rt_gas = (3*Vt_gas/(4*pi))^(1/3);
            t_gas = Pi_gas*Rt_gas/(2*sigma_gas);
            mt_gas = rho_tank_gas*(4*pi/3)*((Rt_gas+t_gas)^3 - Rt_gas^3);
%         else 
%             pressurant_tank.shape = 'cylindrical';
%             % Cylindrical with caps case  
%             h_gas = gas_tank.height;        % Different values of the height of cylinder tank, only for this preliminary phase
% %              N = Vt_real/Vt_gas;        % Number of tanks required (higher or lower than one)
% %              ht_gas = (Vt_gas - (4*pi/3)*Rt_real^3)/(pi*Rt_real^2);         % Pressurant tank heigth [m]
%             f = @(x) Vt_gas - (4*pi/3)*x^3 - pi*h_gas*x^2;
%             Rt_gas = fzero(@(x) f, 0.3*h_gas);
%             t_gas = Pi_gas*Rt_gas/(sigma_gas);
%             mt_gas = rho_tank_gas*((4*pi/3)*((Rt_gas+t_gas)^3 - Rt_gas^3) +...
%                 pi*h_gas*(Rt_gas+t_gas)^2 - Rt_gas^2);
%         end

    case 'Pressure regulated'
        % PROPELLANT
        m_fu = (1/(of+1))*mp;                               % fuel mass
        m_ox = (of/(of+1))*mp;                              % oxidizer mass
        rho_mean = mp/(m_fu/rho_fu + m_ox/rho_ox);          % propellant mean density
        % Margins for ullage
        V_p = (m_fu*1.03+m_ox*1.03)/rho_mean;

%        if V_p < 0.3                                     % Average volume limit assumed from Ariane datasheet
            propellant_tank.shape = 'spherical';
            % Spherical
            Rt_prop = (3*V_p/(4*pi))^(1/3);
            t_prop = Pi_p*Rt_prop/(2*sigma_p);
            mt_prop = rho_tank_p*(4*pi/3)*((Rt_prop+t_prop)^3 - Rt_prop^3);
%        else 
%             propellant_tank.shape = 'cylindrical';
%             % Cylindrical with caps
%             %h_prop = prop_tank.height;        % Different values of the height of cylinder tank, only for this preliminary phase
% 
%             f = @(x) V_p - (4*pi/3)*x^3 - pi*h_prop*x^2;
%             Rt_prop = fzero(@(x) f, 0.3*h_prop);
%             t_prop = Pi_p*Rt_prop/(sigma_p);
%             mt_prop = rho_tank_p*((4*pi/3)*((Rt_prop+t_prop)^3 - Rt_prop^3) +...
%                 pi*h_prop*(Rt_prop+t_prop)^2 - Rt_prop^2);
%         end

        % PRESSURANT
        % Conservation law, adiabatic
        Pi_gas = 10*Pi_p;                                   % Typical assumption, fix the propellant pressure in tanks considering line losses
        m_gas = (Pi_p*V_p*gamma_gas)/(R_gas*T_tank*(1-1/10));      % mass of gas needed
        m_gas = m_gas*1.2;                                         % Margin MAR-MAS-090
        Vt_gas = m_gas*R_gas*T_tank/Pi_gas;                 % Equation of state

%         if Vt_gas < 0.3                                     % Average volume limit assumed from Ariane datasheet
            pressurant_tank.shape = 'spherical';
        % Spherical case
        Rt_gas = (3*Vt_gas/(4*pi))^(1/3);
        t_gas = Pi_gas*Rt_gas/(2*sigma_gas);
        mt_gas = rho_tank_gas*(4*pi/3)*((Rt_gas+t_gas)^3 - Rt_gas^3);
%         else 
%             pressurant_tank.shape = 'cylindrical';
%         % Cylindrical with caps case
% h_gas = gas_tank.height;        % Different values of the height of cylinder tank, only for this preliminary phase% 
%         f = @(x) Vt_gas - (4*pi/3)*x^3 - pi*h_gas*x^2;
%         Rt_gas = fzero(@(x) f, 0.3*h_gas);
%         t_gas = Pi_gas*Rt_gas/(sigma_gas);
%         mt_gas = rho_tank_gas*((4*pi/3)*((Rt_gas+t_gas)^3 - Rt_gas^3) +...
%             pi*h_gas*(Rt_gas+t_gas)^2 - Rt_gas^2);
%     end

end

% Put results in propellant / pressurant tanks structures
propellant_tank.mass = mt_prop; propellant_tank.volume = V_p; 
propellant_tank.radius = Rt_prop;

pressurant_tank.mass = mt_gas; pressurant_tank.volume = Vt_gas; 
pressurant_tank.radius = Rt_gas;

end


