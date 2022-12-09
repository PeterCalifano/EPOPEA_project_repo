function [balances] = HeatBalance_Lander_new(T, R, Q_ext , Q_diss, sigma)

T1 = T(1);
T2 = T(2);
T3 = T(3);
T4 = T(4);

% Balance 1
Q_10 = R.R_10*sigma*(T1^4 - 0);
Q_12 = R.R_21*sigma*(T1^4 - T2^4);
Q_13 = R.R_31*sigma*(T1^4 - T3^4);
Q_14 = R.R_41*sigma*(T1^4 - T4^4);

balances(1) = Q_10 + Q_12 + Q_13 + Q_14 - Q_ext(1) - Q_diss;

% Balance 2
Q_21 = - Q_12;
Q_20 = R.R_20*sigma*(T2^4 - 0);

balances(2) = Q_20 + Q_21 - Q_ext(2);

% Balance 3
Q_31 = - Q_13;
Q_30 = R.R_30*sigma*(T3^4 - 0);

balances(3) = Q_30 + Q_31 - Q_ext(3); 

% Balance 4
Q_41 = - Q_14;
Q_40 = R.R_40*sigma*(T4^4 - 0);

balances(4) = Q_40 + Q_41 - Q_ext(4);


end
