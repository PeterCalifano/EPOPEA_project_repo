function [] = biprop_sizing()

% Add margins according to CDF standards 
V_p = V_p*1.1;                                % Margin MAR-CP-010
% Propellant characteristics
rho_fu = propellant.rho_fu;
rho_ox = propellant.rho_ox;
% Pressurant characteristics
R_gas = pressurant.R;
gamma_gas = pressurant.gamma;
% Pressurant tanks initial properties
Pi_gas = 22*1e+5;               % from datasheet: PROPELLANT TANK FOR SPINSTABILIZED SPACECRAFT T 11/0
% Pressurant tank characteristics
sigma_gas = gas_tank.sigma;
rho_tank_gas = gas_tank.rho;
Rt_real = gas_tank.radius;
Vt_real = gas_tank.volume;

% Assumptions
T_tank = 293;                   % Ambient temperature [K]
switch pressure.mode
    case 'Blowdown'
        m_fu = (1/(of+1))*mp;                               % fuel mass
        m_ox = (of/(of+1))*mp;                              % oxidizer mass
        rho_mean = mp/(m_fu/rho_fu + m_ox/rho_ox);          % propellant mean density
        V_p = mp/rho_mean;
        B = pressure.B;                                     % Chosen blowdown ratio
        Vi_gas = V_p/(B-1);                                 % initial volume of gas needed
        m_gas = (Pi_gas*Vi_gas)/(R_gas*T_tank);
        m_gas = m_gas*1.2;                                  % Margin MAR-MAS-090
        Vt_gas = Vi_gas;                                    % Gas tank volume
        if Vt_gas < 0.3                                     % Average volume limit taken from Ariane datasheet
            % Spherical case
            Rt_gas = (3*Vt_gas/(4*pi))^(1/3);
            t_gas = Pi_gas*Rt_gas/(2*sigma_gas);
            mt_gas = rho_tank_gas*(4*pi/3)*((Rt_gas+t_gas)^3 - Rt_gas^3);
        else 
            % Cylindrical with caps case                    
             N = Vt_real/Vt_gas;        % Number of tanks required (higher or lower than one)
             ht_gas = (Vt_gas - (4*pi/3)*Rt_real^3)/(pi*Rt_real^2);         % Pressurant tank heigth [m]
             t_gas = Pi_gas*Rt_real/(sigma_gas);
        end

    case 'Pressure regulated'


end


end


