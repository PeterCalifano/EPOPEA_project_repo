%% Theoretical RTG sizing
clear ; close all ; clc ;

P_electric_required_watt = 100 ;
t_years = 20 ;
conversion_efficiency = 0.06 ;
try
    fuel_data_Pu = load('Fuel_data_plutonium.mat') ;
    fuel_data_Pu = fuel_data_Pu.Fuel_data ;
    fuel_data_Am = load('Fuel_data_americium.mat') ;
    fuel_data_Am = fuel_data_Am.Fuel_data ;
catch
    error('Remember to add complete EPOPEA repository to path')
end


[ M_isotope_Pu, P_totalThermal_Pu, P_dissipatedThermal_Pu ] = RTGSizing_theoretical( P_electric_required_watt, t_years, conversion_efficiency, fuel_data_Pu )


[ M_isotope_Am, P_totalThermal_Am, P_dissipatedThermal_Am ] = RTGSizing_theoretical( P_electric_required_watt, t_years, conversion_efficiency, fuel_data_Am )


%% Real RTG sizing
clear ; close all ; clc ;

try
    RTG_data_MMRTG = load('RTG_data_MMRTG.mat') ;
    RTG_data_MMRTG = RTG_data_MMRTG.RTG_data ;
    RTG_data_GPHSRTG = load('RTG_data_GPHSRTG.mat') ;
    RTG_data_GPHSRTG = RTG_data_GPHSRTG.RTG_data ;
catch
    error('Remember to add complete EPOPEA repository to path')
end

P_required_EoL = 105 ;
t_EoL_years = 10 ;

[ Number_RTGs_MM, M_tot_RTGs_MM, P_dissipatedThermalTotal_MM, NuclearFuelMassTOTAL_MM ] = RTGSizing_real( P_required_EoL, t_EoL_years, RTG_data_MMRTG )


[ Number_RTGs_GP, M_tot_RTGs_GP, P_dissipatedThermalTotal_GP, NuclearFuelMassTOTAL_GP ] = RTGSizing_real( P_required_EoL, t_EoL_years, RTG_data_GPHSRTG )




