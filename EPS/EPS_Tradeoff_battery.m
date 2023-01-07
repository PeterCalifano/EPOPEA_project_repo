%% Saft VES16 cells -> Big battery (11s16p)
clearvars ; close all ; clc ;
load("Batt_data_LiIon_Saft_VES16_singleBatt.mat") ;
load("Saft_VES16_LiIonCell") ;

% Set inputs
powerRegulationMethod = 'DET' ;
Batt_capacityMargin_percent = 1 ;
PowerDefect_watt = 37.87 ;
PowerDefectDuration_hours = 0.5 ; % [h]
[ TheoreticalBattCapacity_Wh_withmargin, ~, ~, ~, ~ ] = BatterySizing_withRTG( PowerDefect_watt, PowerDefectDuration_hours, Batt_data, Batt_capacityMargin_percent, powerRegulationMethod ) ;

% REFINED ANALYSIS -> Saft cell VES16 (Li-Ion)
V_required = 28 ; % [V] voltage required to power the system
C_required = TheoreticalBattCapacity_Wh_withmargin ; % [Wh] required capacity, coming from the preliminary analysis
[nSeries, nParallel, V_real, C_real] = refinedBatterySizing_withRTG(Cell_data, C_required, V_required) ;

% Compute volume:
D_cell = 0.033 ; % diameter [m]
H_cell = 0.06 ; % height [m]

Tot_volume = D_cell^2*pi/4*H_cell*(nSeries*nParallel) * 1e3 ; % [dm^3]

% Compute mass:
M_cell = 0.155 ; % mass [kg]
Tot_mass = M_cell*(nSeries*nParallel) ;

fprintf( 'Number of batteries (no redundancy): %d\n', Batt_data.NumberOfBatteries ) ;
fprintf( 'Secondary battery voltage [V]: %.2f\n', V_real ) ;
fprintf( 'Secondary battery capacity [Wh]: %.2f\n', C_real ) ;
fprintf( 'Secondary battery volume [dm^3]: %.2f\n', Tot_volume ) ;
fprintf( 'Secondary battery mass [kg]: %.2f\n', Tot_mass ) ;
fprintf( '-------------------------------------------\n' ) ;

Redundancy = 3 ; % Triple redundancy considered on these batteries, due to length of mission

fprintf( 'Number of batteries (with redundancy): %d\n', Batt_data.NumberOfBatteries * Redundancy ) ;
fprintf( 'Total secondary battery capacity [Wh]: %.2f\n', C_real * Redundancy ) ;
fprintf( 'Total secondary battery volume [dm^3]: %.2f\n', Tot_volume * Redundancy ) ;
fprintf( 'Total secondary battery mass [kg]: %.2f\n', Tot_mass * Redundancy ) ;
fprintf( '-------------------------------------------\n' ) ;

%% PRIMARY BATTERY DESIGN Li-SOCl battery
clearvars ; close all ; clc ;

Batt_data.NumberOfBatteries = 3 ;
Batt_data.cycles = 1 ;
Batt_data.DOD = 1 ;
Batt_data.Em = 480 ; % Wh/Kg
Batt_data.Ev = 1024 ; % Wh/dm^3

% Set inputs
powerRegulationMethod = 'DET' ;
Batt_capacityMargin_percent = 0.20 ; % for 15 y (time before landing)
PowerDefect_watt = 535.62 ;
PowerDefectDuration_hours = 1 ; % [h]
[ TheoreticalBattCapacity_Wh_withmargin, TheoreticalBattCapacity_Wh_nomargin, TheoreticalBattMass_kg_withmargin, TheoreticalBattVolume_dm3_withmargin, RequiredEnergy_Wh_nomargin ] = BatterySizing_withRTG( PowerDefect_watt, PowerDefectDuration_hours, Batt_data, Batt_capacityMargin_percent, powerRegulationMethod ) ;

% REFINED ANALYSIS -> PRIMARY BATTERY Li-SOCl2, LS 14250 from Saft
Cell_data.C = 1.2 ;
Cell_data.V = 3.6 ;

V_required = 28 ; % [V] voltage required to power the system
C_required = TheoreticalBattCapacity_Wh_withmargin ; % [Wh] required capacity, coming from the preliminary analysis
[nSeries, nParallel, V_real, C_real] = refinedBatterySizing_withRTG(Cell_data, C_required, V_required) ;

% Compute volume:
D_cell = 0.01462 ; % diameter [m]
H_cell = 0.02513 ; % height [m]

Tot_volume = D_cell^2*pi/4*H_cell*(nSeries*nParallel) * 1e3 ;

% Compute mass:
M_cell = 0.009 ; % mass [kg]
Tot_mass = M_cell*(nSeries*nParallel) ;

Redundancy = 3 ; % Triple redundancy considered on these batteries, due to length of mission

fprintf( 'Number of batteries (no redundancy): %d\n', Batt_data.NumberOfBatteries ) ;
fprintf( 'Primary battery voltage [V]: %.2f\n', V_real ) ;
fprintf( 'Primary battery capacity [Wh]: %.2f\n', C_real ) ;
fprintf( 'Primary battery volume [dm^3]: %.2f\n', Tot_volume ) ;
fprintf( 'Primary battery mass [kg]: %.2f\n', Tot_mass ) ;
fprintf( '-------------------------------------------\n' ) ;

Redundancy = 3 ; % Triple redundancy considered on these batteries, due to length of mission

fprintf( 'Number of batteries (with redundancy): %d\n', Batt_data.NumberOfBatteries * Redundancy ) ;
fprintf( 'Total primary battery capacity [Wh]: %.2f\n', C_real * Redundancy ) ;
fprintf( 'Total primary battery volume [dm^3]: %.2f\n', Tot_volume * Redundancy ) ;
fprintf( 'Total primary battery mass [kg]: %.2f\n', Tot_mass * Redundancy ) ;
fprintf( '-------------------------------------------\n' ) ;
