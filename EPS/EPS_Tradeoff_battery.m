clearvars ; close all ; clc ;
set( 0, 'defaultlegendinterpreter', 'latex' ) ;
set( 0, 'defaulttextinterpreter', 'latex' ) ;

color = { [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840] } ;

%% NOTE: (by Davide)
% First sections are the preliminary analysis with some data taken from
% existing batteries (+ the general case using Lavagna's slides).
% Last two sections refine the sizing using the specific Li-Ion cell of
% those same batteries, so that we can check the consistency of the overall
% Voltage and Capacity [Wh]. All data are taken from Saft and Eaglepicher
% (two producers of batteries for space applications). Datasheets can
% be found in the One Drive PM3/EPS folder.

%% General battery (Li-Ion)
clearvars ; close all ; clc ;
load("Batt_data_LiIon_generic_singleBatt.mat")

% Set inputs
powerRegulationMethod = 'DET';
Batt_capacityMargin_percent = 0.13; 
PowerDefect_watt = 77.3;
PowerDefectDuration_hours = 5690/3600; % [h]
[ TheoreticalBattCapacity_Wh_withmargin, TheoreticalBattMass_kg_withmargin, TheoreticalBattVolume_dm3_withmargin, RequiredEnergy_Wh_nomargin ] = BatterySizing_withRTG( PowerDefect_watt, PowerDefectDuration_hours, Batt_data, Batt_capacityMargin_percent, powerRegulationMethod )

%% Eaglepicher battery
clearvars ; close all ; clc ;
load("Batt_data_LiIon_EP_LP33165_singleBatt.mat")

% Set inputs
powerRegulationMethod = 'DET';
Batt_capacityMargin_percent = 0.13;
PowerDefect_watt = 77.3;
PowerDefectDuration_hours = 5690/3600; % [h]
[ TheoreticalBattCapacity_Wh_withmargin, TheoreticalBattMass_kg_withmargin, TheoreticalBattVolume_dm3_withmargin, RequiredEnergy_Wh_nomargin ] = BatterySizing_withRTG( PowerDefect_watt, PowerDefectDuration_hours, Batt_data, Batt_capacityMargin_percent, powerRegulationMethod )

%% Saft VES16 cells -> Big battery (11s16p)
clearvars ; close all ; clc ;
load("Batt_data_LiIon_Saft_VES16_singleBatt.mat")

% Set inputs
powerRegulationMethod = 'DET';
Batt_capacityMargin_percent = 0.13; 
PowerDefect_watt = 77.3;
PowerDefectDuration_hours = 5690/3600; % [h]
[ TheoreticalBattCapacity_Wh_withmargin, TheoreticalBattMass_kg_withmargin, TheoreticalBattVolume_dm3_withmargin, RequiredEnergy_Wh_nomargin ] = BatterySizing_withRTG( PowerDefect_watt, PowerDefectDuration_hours, Batt_data, Batt_capacityMargin_percent, powerRegulationMethod )

%% Saft VES16 cells -> Small battery (8s4p)
clearvars ; close all ; clc ;
load("Batt_data_LiIon_Saft_VES16_small_singleBatt.mat")

% Set inputs
powerRegulationMethod = 'DET';
Batt_capacityMargin_percent = 0.13; 
PowerDefect_watt = 77.3;
PowerDefectDuration_hours = 5690/3600; % [h]
[ TheoreticalBattCapacity_Wh_withmargin, TheoreticalBattMass_kg_withmargin, TheoreticalBattVolume_dm3_withmargin, RequiredEnergy_Wh_nomargin ] = BatterySizing_withRTG( PowerDefect_watt, PowerDefectDuration_hours, Batt_data, Batt_capacityMargin_percent, powerRegulationMethod )

%% REFINED ANALYSIS -> EaglePicher cells
% load("EP_60Ah_LiIonCell");
% V_required = 28; % [V] voltage required to power the system
% C_required = TheoreticalBattCapacity_Wh_withmargin; % [Wh] required capacity, coming from the preliminary analysis
% [nSeries, nParallel, V_real, C_real] = refinedBatterySizing_withRTG(Cell_data, C_required, V_required)


%% REFINED ANALYSIS -> Saft cell VES16 (Li-Ion)
load("Saft_VES16_LiIonCell") 
V_required = 28; % [V] voltage required to power the system
C_required = TheoreticalBattCapacity_Wh_withmargin; % [Wh] required capacity, coming from the preliminary analysis
[nSeries, nParallel, V_real, C_real] = refinedBatterySizing_withRTG(Cell_data, C_required, V_required)