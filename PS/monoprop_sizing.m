function [] = monoprop_sizing(mp,propellant,pressurant,tank)

% Add margins according to CDF standards 
mp = mp*1.13;                                 % Margin MAR-CP-010 + 3% "unusable" propellant
% Propellant characteristics
rho_p = propellant.rho;
% Tank properties
sigma = tank.sigma;
rho_tank = tank.rho;
h_t = tank.height;        % Different values of the height of cylinder tank, only for this preliminary phase
% Pressurant characteristics
R_gas = pressurant.R; gamma_gas = pressurant.gamma;
% Pressurant tanks initial properties
Pi_gas = 250*1e+5;                            % Change this value! Its an average from datasheets            

% Assumptions
T_tank = 293;                   % Ambient temperature [K]

% BLOWDOWN (propellant and pressurant in the same tank/s)
V_p = mp/rho_p;
B = pressurant.B;                                   % Chosen blowdown ratio
Vt_gas = V_p/(B-1);                                 % initial volume of gas needed
Vt = 1.01*Vt_gas;                                   % Add bladder mass
m_gas = (Pi_gas*Vt_gas)/(R_gas*T_tank);             % Gas initial and tank volume
m_gas = m_gas*1.2;                                  % Margin MAR-MAS-090

% SINGLE TANK SIZING
if Vt < 0.3                                     % Average volume limit assumed from Ariane datasheet
    % Spherical case
    Rt = (3*Vt/(4*pi))^(1/3);
    t = Pi_gas*Rt/(2*sigma);
    mt = rho_tank*(4*pi/3)*((Rt+t)^3 - Rt^3);
else 
    f = @(x) Vt - (4*pi/3)*x^3 - pi*h_t*x^2;
    Rt = fzero(@(x) f, 0.3*h_t);
    t = Pi_gas*Rt/(sigma);
    mt = rho_tank*((4*pi/3)*((Rt+t)^3 - Rt^3) +...
        pi*h_t*(Rt+t)^2 - Rt^2);
end
end

