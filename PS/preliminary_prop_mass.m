function mp = preliminary_prop_mass(dV,m0,Is)
% Preliminary propellant mass sizing
% Computes total propellant mass required for the total dV budget, having
% preliminary estimated the total dry mass at launch
% INPUTS
%       dV = total dV computed (including margins, according to CDF margins
%       standard) [m/s]
%       m0[1x1] = total dry mass at launch (including margins, according to
%       CDF margins standard) [kg]
%
% OUTPUTS
%       mp[1x1] = total propellant mass required for the entire sequence of
%       maneuvres [kg] (Including 2% margin for propellant residuals)
%       m_fu [1x1] = fuel mass [kg] according to selected o/f ratio
%       m_ox [1x1] = oxidizer mass [kg] according to selected o/f ratio
%

% Constants 
g0 = 9.81;                                          % [m/s^2]                           

% Apply Tsiolkowsky equation
mf = m0/exp(dV/(Is*g0));
mp = m0 - mf;
% apply margins 
% from CDF: "the propellant mass shall include 2% for propellant
% residuals."
mp = mp*1.02;

end
