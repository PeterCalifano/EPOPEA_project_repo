clearvars ; close all ; clc ;

% Script generates structure containing solar array data, to loa

SA_data.Id = 0.77 ; % Inherent degradation, area effectively used by solar cells
SA_data.eta_SA = 0.3 ; % Solar panel conversion efficiency
SA_data.yearlyDegradation_percent = 0.0375 ; % Yearly degradation percentage of solar panel efficiency
VarName = 'SA_data_tripleJunction_generic' ;
save( VarName, 'SA_data' ) ;
