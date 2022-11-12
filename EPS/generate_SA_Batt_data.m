%% Generates structure containing solar array data

clearvars ; close all ; clc ;


SA_data.Id = 0.77 ; % Inherent degradation, area effectively used by solar cells
SA_data.eta_SA = 0.3 ; % Solar panel conversion efficiency
SA_data.yearlyDegradation_percent = 0.0375 ; % Yearly degradation percentage of solar panel efficiency
VarName = 'SA_data_tripleJunction_generic' ;
save( VarName, SA_data ) ;

%% Generates structure containing battery data

clearvars ; close all ; clc ;

Batt_data.NumberOfBatteries = 1 ; % Number of batteries
Batt_data.DOD = 0.44 ; % Depth of discharge [%]
Batt_data.cycles = 20000 ; % Battery cycles that can be obtained with the selected DOD
Batt_data.Em = 150 ; % Mass-specific energy [Wh/kg]
Batt_data.Ev = 250 ; % Volume-specific energy [Wh/dm^3]
VarName = 'Batt_data_LiIon_generic_singleBatt' ;
save( VarName, 'Batt_data' ) ;
