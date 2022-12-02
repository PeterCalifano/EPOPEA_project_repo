clear ; close all ; clc ;

clearvars ; close all ; clc ;
set( 0, 'defaultlegendinterpreter', 'latex' ) ;
set( 0, 'defaulttextinterpreter', 'latex' ) ;
color = { [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840] } ;

% Load RTG data
try
    fuel_data_Pu = load('Fuel_data_plutonium.mat') ;
    fuel_data_Pu = fuel_data_Pu.Fuel_data ;
    fuel_data_Am = load('Fuel_data_americium.mat') ;
    fuel_data_Am = fuel_data_Am.Fuel_data ;
    RTG_data_MMRTG = load('RTG_data_MMRTG.mat') ;
    RTG_data_MMRTG = RTG_data_MMRTG.RTG_data ;
    RTG_data_nextGenRTG16 = load('RTG_data_nextGenRTG16.mat') ;
    RTG_data_nextGenRTG16 = RTG_data_nextGenRTG16.RTG_data ;
    RTG_data_nextGenRTG14 = load('RTG_data_nextGenRTG14.mat') ;
    RTG_data_nextGenRTG14 = RTG_data_nextGenRTG14.RTG_data ;
    RTG_data_nextGenRTG12 = load('RTG_data_nextGenRTG12.mat') ;
    RTG_data_nextGenRTG12 = RTG_data_nextGenRTG12.RTG_data ;
    RTG_data_nextGenRTG10 = load('RTG_data_nextGenRTG10.mat') ;
    RTG_data_nextGenRTG10 = RTG_data_nextGenRTG10.RTG_data ;
    RTG_data_nextGenRTG8 = load('RTG_data_nextGenRTG8.mat') ;
    RTG_data_nextGenRTG8 = RTG_data_nextGenRTG8.RTG_data ;
    RTG_data_nextGenRTG6 = load('RTG_data_nextGenRTG6.mat') ;
    RTG_data_nextGenRTG6 = RTG_data_nextGenRTG6.RTG_data ;
    RTG_data_nextGenRTG4 = load('RTG_data_nextGenRTG4.mat') ;
    RTG_data_nextGenRTG4 = RTG_data_nextGenRTG4.RTG_data ;
    RTG_data_GPHSRTG = load('RTG_data_GPHSRTG.mat') ;
    RTG_data_GPHSRTG = RTG_data_GPHSRTG.RTG_data ;
    RTG_data_ASRG = load('RTG_data_ASRG.mat') ;
    RTG_data_ASRG = RTG_data_ASRG.RTG_data ;
catch
    error('Remember to add complete EPOPEA repository to path')
end

% ONLY CHANGE THINGS IN THIS BOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose what architecture to size
% Architecture_choice = 'NSO';
% Architecture_choice = 'SL';
Architecture_choice = 'NSO';

RTGdata1 = RTG_data_nextGenRTG10 ; % This RTG type can be tuned since it is modular, do this once power is known.
RTGdata2 = RTG_data_nextGenRTG8 ; % This RTG type can be tuned since it is modular, do this once power is known.

Number_RTGs1 = 2 ;
Number_RTGs2 = 1 ;

t_EoL_years = 20 ;

% ------------------------------ DO NOT CHANGE BELOW HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute EoL plutonium conditions for each RTG type
[ P_EoL_electric_unit1, P_EoL_electric_unit2 ] = Refined_RTGSizing_real( t_EoL_years, RTGdata1, RTGdata2 ) ;

%% 

% Power required
switch Architecture_choice

    case 'NSO'

        P_req = 400 ;

    case 'SL'

        P_req = 750 ;

    case 'SO'

        P_req = 850 ;

    otherwise

        error('Incorrect architecture entered') ;

end

fprintf( '\nPower required EoL: %.2f', P_req ) ;

[ P_EoL_electric_total, M_tot_RTGs, P_dissipatedThermalTotal, NuclearFuelMassTOTAL ] = Refined_RTGSizing_real_AUXILIARYDATA( t_EoL_years, Number_RTGs1, Number_RTGs2, RTGdata1, RTGdata2 ) ;

fprintf( '\nPower produced EoL: %.2f', P_EoL_electric_total ) ;










