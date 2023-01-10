function [balances] = HeatBalance_Orbiter(T, R,C, Q_ext , Q_diss, Clamped)

T1 = T(1);
T2 = T(2);
T3 = T(3);
T4 = T(4);
T5 = T(5);
T6 = T(6);
Tant = T(7);
Trad = T(8);
T6ext = T(9);
T3ext = T(10);
T2ext = T(11);
T1ext = T(12);
T4ext = T(13);
T5ext = T(14);
T15 = T(15);

% Balance 1 
Q_12 = R.R_12*(T1^4 - T2^4) + C.C_12*(T1 - T2);
Q_13 = R.R_13*(T1^4 - T3^4) + C.C_13*(T1 - T3);
Q_14 = R.R_14*(T1^4 - T4^4) + C.C_14*(T1 - T4);
Q_15 = R.R_15*(T1^4 - T5^4) + C.C_15*(T1 - T5);
Q_16 = R.R_16*(T1^4 - T6^4) + C.C_16*(T1 - T6);
Q_115 = R.R_115*(T1^4 - T15^4) + C.C_115*(T1 - T15);
Q_1ant = R.R_1ant*(T1^4 - Tant^4) + C.C_1ant*(T1 - Tant);
Q_1rad = R.R_1rad*(T1^4 - Trad^4) + C.C_1rad*(T1 - Trad);
Q_10 = R.R_10*(T1ext^4 - 0);
T1_hat = (T1ext + T1)/2;
R1_cond = C.C_1int1ext/4/T1_hat^3;
R1_tot = (1/R.R_1int1ext + 1/R1_cond )^(-1);
Q_1int1ext = R1_tot*(T1^4 - T1ext^4);

balances(1) =  + Q_12 + Q_13 + Q_14  + Q_15  + Q_16 + Q_115 + ...
     + Q_1ant+  + Q_1rad + Q_1int1ext - Q_diss(1);
balances(12) = Q_10  - Q_ext(1) - Q_1int1ext;

% Balance 2
Q_21 = - Q_12;
Q_23 = R.R_23*(T2^4 - T3^4) + C.C_23*(T2 - T3);
Q_24 = R.R_24*(T2^4 - T4^4) + C.C_24*(T2 - T4);
Q_25 = R.R_25*(T2^4 - T5^4) + C.C_25*(T2 - T5);
Q_26 = R.R_26*(T2^4 - T6^4) + C.C_26*(T2 - T6);
Q_2ant = 0;
Q_2rad = R.R_rad2*(T2^4 - Trad^4) + C.C_2rad*(T2 - Trad); 
Q_20 = R.R_20*(T2ext^4 - 0);
T2_hat = (T2ext + T2)/2;
R2_cond = C.C_2int2ext/4/T2_hat^3;
R2_tot = (1/R.R_2int2ext + 1/R2_cond )^(-1);
Q_2int2ext = R2_tot*(T2^4 - T2ext^4);
Q_215 = R.R_215*(T2^4 - T15^4) + C.C_215*(T2 - T15);
balances(11) = Q_20  - Q_ext(2) - Q_2int2ext;

balances(2) = + Q_21 + Q_23 + Q_24  + Q_25  + Q_26 +Q_215+...
     + Q_2ant+  + Q_2rad + Q_2int2ext - Q_diss(2);

% Balance 3
Q_31 = - Q_13;
Q_32 = - Q_23;
Q_34 = R.R_24*(T3^4 - T4^4) + C.C_34*(T3 - T4);
Q_35 = R.R_25*(T3^4 - T5^4) + C.C_35*(T3 - T5);
Q_36 = R.R_26*(T3^4 - T6^4) + C.C_36*(T3 - T6);
Q_315 = R.R_315*(T3^4 - T15^4) + C.C_315*(T3 - T15);
Q_3ant = 0;
Q_3rad = R.R_rad3*(T3^4 - Trad^4) + C.C_3rad*(T3 - Trad); 
Q_30 = R.R_30*(T3ext^4 - 0);
T3_hat = (T3ext + T3)/2;
R3_cond = C.C_3int3ext/4/T3_hat^3;
R3_tot = (1/R.R_3int3ext + 1/R3_cond )^(-1);
Q_3int3ext = R3_tot*(T3^4 - T3ext^4);

balances(3) = + Q_31 + Q_32 + Q_34  + Q_35  + Q_36 + Q_315 + ...
            + Q_3ant+  + Q_3rad + Q_3int3ext - Q_diss(3);
balances(10) = Q_30  - Q_ext(3) - Q_3int3ext;


% Balance 4
Q_41 = - Q_14;
Q_42 = - Q_24;
Q_43 = - Q_34;
Q_45 = R.R_45*(T4^4 - T5^4) + C.C_45*(T4 - T5);
Q_46 = R.R_46*(T4^4 - T6^4) + C.C_46*(T4 - T6);
Q_415 = R.R_415*(T4^4 - T15^4) + C.C_415*(T4- T15);
Q_4ant = 0;
Q_4rad = R.R_rad4*(T4^4 - Trad^4)  + C.C_4rad*(T4 - Trad); 
T4_hat = (T4ext + T4)/2;
R4_cond = C.C_4int4ext/4/T4_hat^3;
R4_tot = (1/R.R_4int4ext + 1/R4_cond )^(-1);
Q_4int4ext = R4_tot*(T4^4 - T4ext^4);

if Clamped ==1
    Q_40 = R.R_40_clamped*(T4ext^4 - 0);
else
    if Clamped == 0
    Q_40 = R.R_40*(T4ext^4 - 0);
    end
end

balances(4) = + Q_41 + Q_42 + Q_43  + Q_45  + Q_46 + Q_415+...
     + Q_4ant+  + Q_4rad + Q_4int4ext - Q_diss(4);
balances(13) = Q_40  - Q_ext(4) - Q_4int4ext;

% Balance 5 
Q_51 = - Q_15;
Q_52 = - Q_25;
Q_53 = - Q_35;
Q_54 = - Q_45;
Q_56 = R.R_56*(T5^4 - T6^4) + C.C_56*(T5 - T6);
Q_515 = R.R_415*(T5^4 - T15^4) + C.C_515*(T5- T15);
Q_5ant = 0;
Q_5rad = R.R_rad5*(T5^4 - Trad^4)  + C.C_5rad*(T5 - Trad);
Q_50 = R.R_50*(T5ext^4 - 0);
T5_hat = (T5ext + T5)/2;
R5_cond = C.C_5int5ext/4/T5_hat^3;
R5_tot = (1/R.R_5int5ext + 1/R5_cond )^(-1);
Q_5int5ext = R5_tot*(T5^4 - T5ext^4);
balances(5) = + Q_51 + Q_52 + Q_53  + Q_54  + Q_56 + Q_515+...
         + Q_5ant+  + Q_5rad - Q_diss(5) + Q_5int5ext;
balances(14) = Q_50  - Q_ext(5) - Q_5int5ext;

% Balance 6 
Q_61 = - Q_16;
Q_62 = - Q_26;
Q_63 = - Q_36;
Q_64 = - Q_46;
Q_65 = - Q_56;
Q_6ant = 0;
Q_615 = R.R_615*(T6^4 - T15^4) + C.C_615*(T6- T15);
Q_6rad = R.R_rad6*(T6^4 - Trad^4)  + C.C_6rad*(T6 - Trad);
Q_6ext0 = R.R_60*(T6ext^4 - 0);
T6_hat = (T6ext + T6)/2;
R6_cond = C.C_6int6ext/4/T6_hat^3;
R6_tot = (1/R.R_6int6ext + 1/R6_cond )^(-1);
Q_6int6ext = R6_tot*(T6^4 - T6ext^4);
balances(9) =  + Q_6ext0 - Q_6int6ext - Q_ext(6);

balances(6) = + Q_61 + Q_62 + Q_63  + Q_64  + Q_65 + Q_615 + ...
     + Q_6ant+  + Q_6rad + Q_6int6ext - Q_diss(6);

% Balance 7 (antenna)
Q_ant1 = - Q_1ant;
Q_ant0 = R.R_ant0*(Tant^4 - 0);

balances(7) = Q_ant0 + Q_ant1- Q_ext(7) - Q_diss(7);

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


% Balance 15

Q_151 = - Q_115;
Q_152 = - Q_215;
Q_153 = - Q_315;
Q_154 = - Q_415;
Q_155 = - Q_515;
Q_156 = - Q_615;
Q_15rad = R.R_rad15*(T15^4 - Trad^4)  + C.C_15rad*(T15 - Trad);
balances(15) = + Q_151 + Q_152 + Q_153  + Q_154  + Q_155 + Q_156 + ...
     + Q_15rad+ - Q_diss(15);

end
