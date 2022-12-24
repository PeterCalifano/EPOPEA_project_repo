%clear; clc; close all
%T_Mimas = 0.9424218;
T_Enceladus = 1.370218;
T_Tethys = 1.887802; 
T_Dione = 2.736915;
T_Rhea = 4.517500;
T_Titan = 15.945421;

T_ET = 1/abs(1/T_Enceladus - 1/T_Tethys);
T_ETD = 1/abs(1/T_ET - 1/T_Dione);
T_ETDR = 1/abs(1/T_ETD - 1/T_Rhea);
T_ETDRT = 1/abs(1/T_ETDR - 1/T_Titan);