function [balances] = HeatBalance_Clamped(T, RO,Rl, CO,Cl, Q_extO ,Q_extl , Q_dissO,Q_dissl, Clamped)

TO = T(1:15);
Tl = T(16:end);

T1o = TO(1);
T2o = TO(2);
T3o = TO(3);
T4o = TO(4);
T5o = TO(5);
T6o = TO(6);
Tant = TO(7);
Trado = TO(8);
T6exto = TO(9);
T3exto = TO(10);
T2exto = TO(11);
T1exto = TO(12);
T4exto = TO(13);
T5exto = TO(14);
T15o = TO(15);

T1l = Tl(1);
T2l = Tl(2);
T3l = Tl(3);
T4l = Tl(4);
T5l = Tl(5);
T6l = Tl(6);
T8l = Tl(8);
Tradl = Tl(7);
T1extl = Tl(10);
T3extl = Tl(9);
T2extl = Tl(11);
T4extl = Tl(12);
T5extl = Tl(13);
T6extl = Tl(14);
% Balance 1 
Q_12 = RO.R_12*(T1o^4 - T2o^4) + CO.C_12*(T1o - T2o);
Q_13 = RO.R_13*(T1o^4 - T3o^4) + CO.C_13*(T1o - T3o);
Q_14 = RO.R_14*(T1o^4 - T4o^4) + CO.C_14*(T1o - T4o);
Q_15 = RO.R_15*(T1o^4 - T5o^4) + CO.C_15*(T1o - T5o);
Q_16 = RO.R_16*(T1o^4 - T6o^4) + CO.C_16*(T1o - T6o);
Q_115 = RO.R_115*(T1o^4 - T15o^4) + CO.C_115*(T1o - T15o);
Q_1ant = RO.R_1ant*(T1o^4 - Tant^4) + CO.C_1ant*(T1o - Tant);
Q_1rad = RO.R_1rad*(T1o^4 - Trado^4) + CO.C_1rad*(T1o - Trado);
Q_10 = RO.R_10*(T1exto^4 - 0);
T1_hat = (T1exto + T1o)/2;
R1_cond = CO.C_1int1ext/4/T1_hat^3;
R1_tot = (1/RO.R_1int1ext + 1/R1_cond )^(-1);
Q_1int1ext = R1_tot*(T1o^4 - T1exto^4);

balanceso(1) =  + Q_12 + Q_13 + Q_14  + Q_15  + Q_16 + Q_115 + ...
     + Q_1ant+  + Q_1rad + Q_1int1ext - Q_dissO(1);
balanceso(12) = Q_10  - Q_extO(1) - Q_1int1ext;

% Balance 2
Q_21 = - Q_12;
Q_23 = RO.R_23*(T2o^4 - T3o^4) + CO.C_23*(T2o - T3o);
Q_24 = RO.R_24*(T2o^4 - T4o^4) + CO.C_24*(T2o - T4o);
Q_25 = RO.R_25*(T2o^4 - T5o^4) + CO.C_25*(T2o - T5o);
Q_26 = RO.R_26*(T2o^4 - T6o^4) + CO.C_26*(T2o - T6o);
Q_2ant = 0;
Q_2rad = RO.R_rad2*(T2o^4 - Trado^4) + CO.C_2rad*(T2o - Trado); 
Q_20 = RO.R_20*(T2exto^4 - 0);
T2_hat = (T2exto + T2o)/2;
R2_cond = CO.C_2int2ext/4/T2_hat^3;
R2_tot = (1/RO.R_2int2ext + 1/R2_cond )^(-1);
Q_2int2ext = R2_tot*(T2o^4 - T2exto^4);
Q_215 = RO.R_215*(T2o^4 - T15o^4) + CO.C_215*(T2o - T15o);
balanceso(11) = Q_20  - Q_extO(2) - Q_2int2ext;

balanceso(2) = + Q_21 + Q_23 + Q_24  + Q_25  + Q_26 +Q_215+...
     + Q_2ant+  + Q_2rad + Q_2int2ext - Q_dissO(2);

% Balance 3
Q_31 = - Q_13;
Q_32 = - Q_23;
Q_34 = RO.R_24*(T3o^4 - T4o^4) + CO.C_34*(T3o - T4o);
Q_35 = RO.R_25*(T3o^4 - T5o^4) + CO.C_35*(T3o - T5o);
Q_36 = RO.R_26*(T3o^4 - T6o^4) + CO.C_36*(T3o - T6o);
Q_315 = RO.R_315*(T3o^4 - T15o^4) + CO.C_315*(T3o - T15o);
Q_3ant = 0;
Q_3rad = RO.R_rad3*(T3o^4 - Trado^4) + CO.C_3rad*(T3o - Trado); 
Q_30 = RO.R_30*(T3exto^4 - 0);
T3_hat = (T3exto + T3o)/2;
R3_cond = CO.C_3int3ext/4/T3_hat^3;
R3_tot = (1/RO.R_3int3ext + 1/R3_cond )^(-1);
Q_3int3ext = R3_tot*(T3o^4 - T3exto^4);

balanceso(3) = + Q_31 + Q_32 + Q_34  + Q_35  + Q_36 + Q_315 + ...
            + Q_3ant+  + Q_3rad + Q_3int3ext - Q_dissO(3);
balanceso(10) = Q_30  - Q_extO(3) - Q_3int3ext;


% Balance 4
Q_41 = - Q_14;
Q_42 = - Q_24;
Q_43 = - Q_34;
Q_45 = RO.R_45*(T4o^4 - T5o^4) + CO.C_45*(T4o - T5o);
Q_46 = RO.R_46*(T4o^4 - T6o^4) + CO.C_46*(T4o - T6o);
Q_415 = RO.R_415*(T4o^4 - T15o^4) + CO.C_415*(T4o- T15o);
Q_4ant = 0;
Q_4rad = RO.R_rad4*(T4o^4 - Trado^4)  + CO.C_4rad*(T4o - Trado); 
T4_hat = (T4exto + T4o)/2;
R4_cond = CO.C_4int4ext/4/T4_hat^3;
R4_tot = (1/RO.R_4int4ext + 1/R4_cond )^(-1);
Q_4int4ext = R4_tot*(T4o^4 - T4exto^4);

if Clamped ==1
    Q_40 = RO.R_40_clamped*(T4exto^4 - 0);
    Q_4o2l = RO.R_4o2l*(T4o^4 - T2l^4) + CO.R_4o2l*(T4o- T2l);
else
    if Clamped == 0
    Q_40 = RO.R_40*(T4exto^4 - 0);
    Q_4o2l = 0;
    end
end

balanceso(4) = + Q_41 + Q_42 + Q_43  + Q_45  + Q_46 + Q_415+...
     + Q_4ant+  + Q_4rad + Q_4int4ext - Q_dissO(4) + Q_4o2l;
balanceso(13) = Q_40  - Q_extO(4) - Q_4int4ext;

% Balance 5 
Q_51 = - Q_15;
Q_52 = - Q_25;
Q_53 = - Q_35;
Q_54 = - Q_45;
Q_56 = RO.R_56*(T5o^4 - T6o^4) + CO.C_56*(T5o - T6o);
Q_515 = RO.R_415*(T5o^4 - T15o^4) + CO.C_515*(T5o- T15o);
Q_5ant = 0;
Q_5rad = RO.R_rad5*(T5o^4 - Trado^4)  + CO.C_5rad*(T5o - Trado);
Q_50 = RO.R_50*(T5exto^4 - 0);
T5_hat = (T5exto + T5o)/2;
R5_cond = CO.C_5int5ext/4/T5_hat^3;
R5_tot = (1/RO.R_5int5ext + 1/R5_cond )^(-1);
Q_5int5ext = R5_tot*(T5o^4 - T5exto^4);
balanceso(5) = + Q_51 + Q_52 + Q_53  + Q_54  + Q_56 + Q_515+...
         + Q_5ant+  + Q_5rad - Q_dissO(5) + Q_5int5ext;
balanceso(14) = Q_50  - Q_extO(5) - Q_5int5ext;

% Balance 6 
Q_61 = - Q_16;
Q_62 = - Q_26;
Q_63 = - Q_36;
Q_64 = - Q_46;
Q_65 = - Q_56;
Q_6ant = 0;
Q_615 = RO.R_615*(T6o^4 - T15o^4) + CO.C_615*(T6o- T15o);
Q_6rad = RO.R_rad6*(T6o^4 - Trado^4)  + CO.C_6rad*(T6o - Trado);
Q_6ext0 = RO.R_60*(T6exto^4 - 0);
T6_hat = (T6exto + T6o)/2;
R6_cond = CO.C_6int6ext/4/T6_hat^3;
R6_tot = (1/RO.R_6int6ext + 1/R6_cond )^(-1);
Q_6int6ext = R6_tot*(T6o^4 - T6exto^4);
balanceso(9) =  + Q_6ext0 - Q_6int6ext - Q_extO(6);

balanceso(6) = + Q_61 + Q_62 + Q_63  + Q_64  + Q_65 + Q_615 + ...
     + Q_6ant+  + Q_6rad + Q_6int6ext - Q_dissO(6);

% Balance 7 (antenna)
Q_ant1 = - Q_1ant;
Q_ant0 = RO.R_ant0*(Tant^4 - 0);

balanceso(7) = Q_ant0 + Q_ant1- Q_extO(7) - Q_dissO(7);

% Balance 8 (rad)
Q_rad1 = - Q_1rad;
Q_rad2 = - Q_2rad;
Q_rad3 = - Q_3rad;
Q_rad4 = - Q_4rad;
Q_rad5 = - Q_5rad;
Q_rad6 = - Q_6rad;
Q_rad0 = RO.R_rad0*(Trado^4 - 0);

balanceso(8) = Q_rad0 + Q_rad1 + Q_rad2 + Q_rad3  + Q_rad4  + Q_rad6 +...
       + Q_rad5- Q_extO(8);


% Balance 15

Q_151 = - Q_115;
Q_152 = - Q_215;
Q_153 = - Q_315;
Q_154 = - Q_415;
Q_155 = - Q_515;
Q_156 = - Q_615;
Q_15rad = RO.R_rad15*(T15o^4 - Trado^4)  + CO.C_15rad*(T15o - Trado);
balanceso(15) = + Q_151 + Q_152 + Q_153  + Q_154  + Q_155 + Q_156 + ...
     + Q_15rad+ - Q_dissO(15);

% Balance 1 
Q_12 = Rl.R_12*(T1l^4 - T2l^4) + Cl.C_12*(T1l - T2l);
Q_13 = Rl.R_13*(T1l^4 - T3l^4) + Cl.C_13*(T1l - T3l);
Q_14 = Rl.R_14*(T1l^4 - T4l^4) + Cl.C_14*(T1l - T4l);
Q_15 = Rl.R_15*(T1l^4 - T5l^4) + Cl.C_15*(T1l - T5l);
Q_16 = Rl.R_16*(T1l^4 - T6l^4) + Cl.C_16*(T1l - T6l);
Q_18 = Rl.R_18*(T1l^4 - T8l^4) + Cl.C_18*(T1l - T8l);
Q_1rad = Rl.R_1rad*(T1l^4 - Tradl^4) + Cl.C_1rad*(T1l - Tradl);
Q_10 = Rl.R_10*(T1extl^4 - 0);
T1_hat = (T1extl + T1l)/2;
R1_cond = Cl.C_1int1ext/4/T1_hat^3;
R1_tot = (1/Rl.R_1int1ext + 1/R1_cond )^(-1);
Q_1int1ext = R1_tot*(T1l^4 - T1extl^4);
balancesl(10) = Q_10  - Q_extl(1) - Q_1int1ext;

balancesl(1) =  + Q_12 + Q_13 + Q_14  + Q_15  + Q_16 +...
       + Q_1rad + Q_18 +  Q_1int1ext - Q_dissl(1);

% Balance 2
Q_21 = - Q_12;
Q_23 = Rl.R_23*(T2l^4 - T3l^4) + Cl.C_23*(T2l - T3l);
Q_24 = Rl.R_24*(T2l^4 - T4l^4) + Cl.C_24*(T2l - T4l);
Q_25 = Rl.R_25*(T2l^4 - T5l^4) + Cl.C_25*(T2l - T5l);
Q_26 = Rl.R_26*(T2l^4 - T6l^4) + Cl.C_26*(T2l - T6l);
Q_28 = Rl.R_28*(T2l^4 - T8l^4) + Cl.C_28*(T2l - T8l);
Q_2rad = Rl.R_rad2*(T2l^4 - Tradl^4) + Cl.C_2rad*(T2l - Tradl); 
T2_hat = (T2extl + T2l)/2;
R2_cond = Cl.C_2int2ext/4/T2_hat^3;
R2_tot = (1/Rl.R_2int2ext + 1/R2_cond )^(-1);
Q_2int2ext = R2_tot*(T2l^4 - T2extl^4);
if Clamped ==1
    Q_20 = 0;
    Q_4l2o = - Q_4o2l;
else
    if Clamped == 0
    Q_20 = Rl.R_20*(T2extl^4 - 0);
    Q_4l2o = 0;
    end
end

balancesl(11) = Q_20  - Q_extl(2) - Q_2int2ext;
balancesl(2) =+ Q_21 + Q_23 + Q_24  + Q_25  + Q_26 +...
      + Q_2rad  + Q_28 - Q_dissl(2) + Q_2int2ext + Q_4l2o;

% Balance 3
Q_31 = - Q_13;
Q_32 = - Q_23;
Q_34 = Rl.R_24*(T3l^4 - T4l^4) + Cl.C_34*(T3l - T4l);
Q_35 = Rl.R_25*(T3l^4 - T5l^4) + Cl.C_35*(T3l - T5l);
Q_36 = Rl.R_26*(T3l^4 - T6l^4) + Cl.C_36*(T3l - T6l);
Q_38 = Rl.R_28*(T3l^4 - T8l^4) + Cl.C_38*(T3l - T8l);
Q_3rad = Rl.R_rad3*(T3l^4 - Tradl^4) + Cl.C_3rad*(T3l - Tradl); 
Q_30 = Rl.R_30*(T3extl^4 - 0);
T3_hat = (T3extl + T3l)/2;
R3_cond = Cl.C_3int3ext/4/T3_hat^3;
R3_tot = (1/Rl.R_3int3ext + 1/R3_cond )^(-1);
Q_3int3ext = R3_tot*(T3l^4 - T3extl^4);

balancesl(9) = Q_30  - Q_extl(3) - Q_3int3ext;

balancesl(3) =  + Q_31 + Q_32 + Q_34  + Q_35  + Q_36 +...
              + Q_3rad + Q_38 + Q_3int3ext - Q_dissl(3);

% Balance 4
Q_41 = - Q_14;
Q_42 = - Q_24;
Q_43 = - Q_34;
Q_45 = Rl.R_45*(T4l^4 - T5l^4) + Cl.C_45*(T4l - T5l);
Q_46 = Rl.R_46*(T4l^4 - T6l^4) + Cl.C_46*(T4l - T6l);
Q_48 = Rl.R_48*(T4l^4 - T8l^4) + Cl.C_28*(T4l - T8l);
Q_4rad = Rl.R_rad4*(T4l^4 - Tradl^4)  + Cl.C_4rad*(T4l - Tradl); 
Q_40 = Rl.R_40*(T4extl^4 - 0);
T4_hat = (T4extl + T4l)/2;
R4_cond = Cl.C_4int4ext/4/T4_hat^3;
R4_tot = (1/Rl.R_4int4ext + 1/R4_cond )^(-1);
Q_4int4ext = R4_tot*(T4l^4 - T4extl^4);
balancesl(4) =  + Q_41 + Q_42 + Q_43  + Q_45  + Q_46 +...
       + Q_4rad+ Q_48 - Q_dissl(4) + Q_4int4ext;
balancesl(12) = Q_40  - Q_extl(4) - Q_4int4ext;
% Balance 5 
Q_51 = - Q_15;
Q_52 = - Q_25;
Q_53 = - Q_35;
Q_54 = - Q_45;
Q_56 = Rl.R_56*(T5l^4 - T6l^4) + Cl.C_56*(T5l - T6l);
Q_58 = Rl.R_58*(T5l^4 - T8l^4) + Cl.C_58*(T5l - T8l);
Q_5rad = Rl.R_rad5*(T5l^4 - Tradl^4)  + Cl.C_5rad*(T5l - Tradl);
Q_50 = Rl.R_50*(T5extl^4 - 0);
T5_hat = (T5extl + T5l)/2;
R5_cond = Cl.C_5int5ext/4/T5_hat^3;
R5_tot = (1/Rl.R_5int5ext + 1/R5_cond )^(-1);
Q_5int5ext = R5_tot*(T5l^4 - T5extl^4);
balancesl(13) = Q_50  - Q_extl(5) - Q_5int5ext;
balancesl(5) = + Q_51 + Q_52 + Q_53  + Q_54  + Q_56 +...
           + Q_5rad   - Q_dissl(5) + Q_58 + Q_5int5ext;

% Balance 6 
Q_61 = - Q_16;
Q_62 = - Q_26;
Q_63 = - Q_36;
Q_64 = - Q_46;
Q_65 = - Q_56;
Q_6rad = Rl.R_rad6*(T6l^4 - Tradl^4)  + Cl.C_6rad*(T6l - Tradl);
Q_60 = Rl.R_60*(T6extl^4 - 0);
Q_68 = Rl.R_28*(T6l^4 - T8l^4) + Cl.C_28*(T6l - T8l);
T6_hat = (T6extl + T6l)/2;
R6_cond = Cl.C_6int6ext/4/T6_hat^3;
R6_tot = (1/Rl.R_6int6ext + 1/R6_cond )^(-1);
Q_6int6ext = R6_tot*(T6l^4 - T6extl^4);
balancesl(14) =  + Q_60 - Q_6int6ext - Q_extl(6);
balancesl(6) = + Q_61 + Q_62 + Q_63  + Q_64  + Q_65 +...
              + Q_6rad  + Q_68 - Q_dissl(6)  + Q_6int6ext;
% Q_6ext0 = R.R_60*(T6ext^4 - 0);
% Q_6int6ext = R.R_6int6ext*(T6^4 - T6ext^4);
% 
% balances(9) =  + Q_6ext0 - Q_6int6ext - Q_ext(6);
% 
% balances(6) = + Q_61 + Q_62 + Q_63  + Q_64  + Q_65 +...
%     + Q_6rad + Q_6int6ext;

% Balance 7 (rad)
Q_rad1 = - Q_1rad;
Q_rad2 = - Q_2rad;
Q_rad3 = - Q_3rad;
Q_rad4 = - Q_4rad;
Q_rad5 = - Q_5rad;
Q_rad6 = - Q_6rad;
Q_rad0 = Rl.R_rad0*(Tradl^4 - 0);
Q_rad8 = Rl.R_rad8*(Tradl^4 - T8l^4) + Cl.C_8rad*(Tradl - T8l);
balancesl(7) = Q_rad0 + Q_rad1 + Q_rad2 + Q_rad3  + Q_rad4  + Q_rad6 +...
       + Q_rad5- Q_extl(7) + Q_rad8 - Q_dissl(7);

% Balance 8 (int node)
Q_81 = - Q_18;
Q_82 = - Q_28;
Q_83 = - Q_38;
Q_84 = - Q_48;
Q_85 = - Q_58;
Q_86 = - Q_68;
Q_8rad = - Q_rad8;
balancesl(8) = + Q_81 + Q_82 + Q_83  + Q_84  + Q_85 +...
              + Q_8rad + Q_86 - Q_dissl(8);

balances = [balanceso, balancesl];

end
