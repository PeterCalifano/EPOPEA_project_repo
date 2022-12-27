function [balances] = HeatBalance_Orbiter(T, R,C, Q_ext , Q_diss, sigma)

T1 = T(1);
T2 = T(2);
T3 = T(3);
T4 = T(4);
T5 = T(5);
T6 = T(6);
Tant = T(7);
Trad = T(8);

% Conductive coupling
% Cond_15 = C.C_15/(4*0.5^3*(T1 + T5)^3);

% Balance 1
Q_12 = R.R_12*(T1^4 - T2^4);
Q_13 = R.R_13*(T1^4 - T3^4);
Q_14 = R.R_14*(T1^4 - T4^4);
Q_15 = R.R_15*(T1^4 - T5^4);
Q_16 = R.R_16*(T1^4 - T6^4);
Q_1ant = R.R_1ant*(T1^4 - Tant^4);
Q_1rad = R.R_1rad*(T1^4 - Trad^4);
Q_10 = R.R_10*(T1^4 - 0);

balances(1) = Q_10 + Q_12 + Q_13 + Q_14  + Q_15  + Q_16 +...
     + Q_1ant+  + Q_1rad- Q_ext(1) - Q_diss;

% Balance 2
Q_21 = - Q_12;
Q_23 = R.R_23*(T2^4 - T3^4);
Q_24 = R.R_24*(T2^4 - T4^4);
Q_25 = R.R_25*(T2^4 - T5^4);
Q_26 = R.R_26*(T2^4 - T6^4);
Q_2ant = 0;
Q_2rad = 0;
Q_20 = R.R_20*(T2^4 - 0);

balances(2) = Q_20 + Q_21 + Q_23 + Q_24  + Q_25  + Q_26 +...
     + Q_2ant+  + Q_2rad- Q_ext(2);

% Balance 3
Q_31 = - Q_13;
Q_32= - Q_23;
Q_34 = R.R_24*(T3^4 - T4^4);
Q_35 = R.R_25*(T3^4 - T5^4);
Q_36 = R.R_26*(T3^4 - T6^4);
Q_3ant = 0;
Q_3rad = R.R_rad3*(T3^4 - Trad^4);
Q_30 = R.R_30*(T3^4 - 0);

balances(3) = Q_30 + Q_31 + Q_32 + Q_34  + Q_35  + Q_36 +...
     + Q_3ant+  + Q_3rad- Q_ext(3);
% Balance 4
Q_41 = - Q_14;
Q_42= - Q_24;
Q_43 = - Q_34;
Q_45 = R.R_45*(T4^4 - T5^4);
Q_46 = R.R_46*(T4^4 - T6^4);
Q_4ant = 0;
Q_4rad = R.R_rad4*(T4^4 - Trad^4);
Q_40 = R.R_40*(T4^4 - 0);

balances(4) = Q_40 + Q_41 + Q_42 + Q_43  + Q_45  + Q_46 +...
     + Q_4ant+  + Q_4rad- Q_ext(4);
% Balance 5 
Q_51 = - Q_15;
Q_52 = - Q_25;
Q_53 = - Q_35;
Q_54 = - Q_45;
Q_56 = R.R_56*(T5^4 - T6^4);
Q_5ant = 0;
Q_5rad = R.R_rad5*(T5^4 - Trad^4);
Q_50 = R.R_50*(T5^4 - 0);

balances(5) = Q_50 + Q_51 + Q_52 + Q_53  + Q_54  + Q_56 +...
     + Q_5ant+  + Q_5rad- Q_ext(5);

% Balance 6 
Q_61 = - Q_16;
Q_62 = - Q_26;
Q_63 = - Q_36;
Q_64 = - Q_46;
Q_65 = - Q_56;
Q_6ant = 0;
Q_6rad = R.R_rad6*(T6^4 - Trad^4);
Q_60 = R.R_60*(T6^4 - 0);

balances(6) = Q_60 + Q_61 + Q_62 + Q_63  + Q_64  + Q_65 +...
     + Q_6ant+  + Q_6rad- Q_ext(6);

% Balance 7 (antenna)
Q_ant1 = - Q_1ant;
Q_ant0 = R.R_ant0*(Tant^4 - 0);

balances(7) = Q_ant0 + Q_ant1- Q_ext(7);

% Balance 8 (rad)
Q_rad1 = - Q_1rad;
Q_rad2 = - Q_2rad;
Q_rad3 = - Q_3rad;
Q_rad4 = - Q_4rad;
Q_rad5 = - Q_5rad;
Q_rad6 = - Q_6rad;
Q_rad0 = R.R_rad0*(Trad^4 - 0);

balances(8) = Q_rad0 + Q_rad1 + Q_rad2 + Q_rad3  + Q_rad4  + Q_rad6 +...
       + Q_rad5- Q_ext(8);

end
