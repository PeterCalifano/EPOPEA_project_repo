clear all, close all, clc
% ----------------------------------- 400 N thruster
R = 2077.3 ;
m_dot = 0.135 ;
rho = 0.178 ; 
m = 11.38 ;
p = 300e5 ;
V = 0.23 ;
T = 25 + 273.15 ;

dp = dpressure( m_dot, rho, m, R, p, V, T )

%% ////////////
function dp = dpressure( m_dot, rho, m, R, p, V, T )

% m_dot = mass flow rate of helium - assumed to be constant and equal to the one required by thruster
% rho = density of helium
% m = mass of helium
% p = initial pressure of helium
% V = initial volume of helium
% T = temperature tank

Q_vol = m_dot / rho ; % [m^3/s] - Volumetric flow (m_dot assumed constant)

PM = 0.004003 ; % [kg/mol] - Helium molecular mass
n = m / PM ; % [mol] 

R = 8.31 ;
dp = ( n * R * T  - Q_vol * p ) / ( V ) ;

end
