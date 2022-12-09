function [balances] = HeatBalance_Lander(T, R, T_E, Q_Sat_T,Q_Sat_L , Q_diss, Q_RTG, sigma,C_LB,C_LT)

T_B = T(1);
T_L = T(2);
T_T = T(3);

Cond_LB = C_LB/(4*sigma*0.5^3*(T_L + T_B)^3);
Cond_LT = C_LT/(4*sigma*0.5^3*(T_L + T_T)^3);

% balance bottom
Q_EB = R.R_EB*sigma*(T_E^4-T_B^4);
Q_LB = (R.R_LB + Cond_LB)*sigma*(T_L^4-T_B^4);
Q_TB = R.R_TB*sigma*(T_T^4-T_B^4);

% balance lateral surfaces
Q_BL = (R.R_BL + Cond_LB) *sigma*(T_B^4-T_L^4);
Q_TL = (R.R_TL + Cond_LT) *sigma*(T_T^4-T_L^4);
Q_LDS = R.R_LDS*sigma*T_L^4;

% balance top surface
Q_BT = R.R_BT*sigma*(T_B^4-T_T^4);
Q_LT = (R.R_LT + Cond_LT)*sigma*(T_L^4-T_T^4);
Q_TDS = R.R_TDS*sigma*T_T^4;

balances(1) = Q_EB+Q_LB+Q_TB;
balances(2) = Q_BL+Q_TL-Q_LDS+Q_Sat_L+Q_diss+Q_RTG;
balances(3) = Q_BT+Q_LT-Q_TDS+Q_Sat_T; 

end
