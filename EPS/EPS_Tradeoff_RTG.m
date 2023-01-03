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

    error('Remember to add complete EPOPEA repository to path') ;

end

% Real RTG data - Used to select which type of RTG to consider
RTG_data = RTG_data_nextGenRTG16 ; % This RTG type can be tuned since it is modular, do this once power is known.
%RTG_data = RTG_data_GPHSRTG ; 
%RTG_data = RTG_data_MMRTG ; 
%RTG_data = RTG_data_ASRG ; 

% Power requirements EoL
%P_req_NSOSL_orb = 375 ;
P_req_lan = 201 ;
P_req_SOSL_orb = 336 ;

% EoL
t_EoL_years_NSOSL_orb = 20 ;
t_EoL_years_NSOSL_lan = 20 ;
t_EoL_years_SOSL_orb = 20 ;
t_EoL_years_SOSL_lan = 20 ;

% %% Compute theoretical RTG sizing for NSOSL (orbiter) architecture
% % --------------------------------------
% 
% P_electric_required_watt_NSOSL_orb = P_req_NSOSL_orb ;
% conversion_efficiency = 0.07 ;
% 
% [ M_isotope_NSOSL_orb_Pu, ~, ~ ] = RTGSizing_theoretical( P_electric_required_watt_NSOSL_orb, t_EoL_years_NSOSL_orb, conversion_efficiency, fuel_data_Pu ) ;
% [ M_isotope_NSOSL_orb_Am, ~, ~ ] = RTGSizing_theoretical( P_electric_required_watt_NSOSL_orb, t_EoL_years_NSOSL_orb, conversion_efficiency, fuel_data_Am ) ;
% 
% 
% %% Compute theoretical RTG sizing for NSOSL (lander) architecture
% % --------------------------------------
% 
% P_electric_required_watt_NSOSL_lan = P_req_lan ;
% conversion_efficiency = 0.07 ;
% 
% [ M_isotope_NSOSL_lan_Pu, ~, ~ ] = RTGSizing_theoretical( P_electric_required_watt_NSOSL_lan, t_EoL_years_NSOSL_lan, conversion_efficiency, fuel_data_Pu ) ;
% [ M_isotope_NSOSL_lan_Am, ~, ~ ] = RTGSizing_theoretical( P_electric_required_watt_NSOSL_lan, t_EoL_years_NSOSL_lan, conversion_efficiency, fuel_data_Am ) ;
% 
% 
% %% Compute theoretical RTG sizing for SO + SL (orbiter) architecture 
% % ------------------------------------------
% 
% P_electric_required_watt_SOSL_orb = P_req_SOSL_orb ;
% conversion_efficiency = 0.07 ;
% 
% [ M_isotope_SOSL_orb_Pu, ~, ~ ] = RTGSizing_theoretical( P_electric_required_watt_SOSL_orb, t_EoL_years_SOSL_orb, conversion_efficiency, fuel_data_Pu ) ;
% [ M_isotope_SOSL_orb_Am, ~, ~ ] = RTGSizing_theoretical( P_electric_required_watt_SOSL_orb, t_EoL_years_SOSL_orb, conversion_efficiency, fuel_data_Am ) ;


% %% Compute real RTG sizing for NSOSL (orbiter) architecture
% % --------------------------------------
% 
% P_required_EoL_NSOSL_orb = P_req_NSOSL_orb ;
% 
% [ Number_RTGs_NSOSL_orb, M_tot_RTGs_NSOSL_orb, P_BOLdissipatedThermalTotal_NSOSL_orb, NuclearFuelMassTOTAL_NSOSL_orb, P_EoL_electric_total_NSOSL_orb ] = RTGSizing_real( P_required_EoL_NSOSL_orb, t_EoL_years_NSOSL_orb, RTG_data ) ;


%% Compute real RTG sizing for NSOSL (lander) architecture
% --------------------------------------

P_required_EoL_NSOSL_lan = P_req_lan ;

[ Number_RTGs_NSOSL_lan, M_tot_RTGs_NSOSL_lan, P_BOLdissipatedThermalTotal_NSOSL_lan, NuclearFuelMassTOTAL_NSOSL_lan, P_EoL_electric_total_NSOSL_lan ] = RTGSizing_real( P_required_EoL_NSOSL_lan, t_EoL_years_NSOSL_lan, RTG_data ) ;


%% Compute real RTG sizing for SO + SL (orbiter) architecture 
% ------------------------------------------

P_required_EoL_SOSL_orb = P_req_SOSL_orb ;

[ Number_RTGs_SOSL_orb, M_tot_RTGs_SOSL_orb, P_BOLdissipatedThermalTotal_SOSL_orb, NuclearFuelMassTOTAL_SOSL_orb, P_EoL_electric_total_SOSL_orb ] = RTGSizing_real( P_required_EoL_SOSL_orb, t_EoL_years_SOSL_orb, RTG_data ) ;

