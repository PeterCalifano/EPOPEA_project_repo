clearvars ; close all ; clc ;


PowerDefect_watt = 1031-660 ;
PowerDefectDuration_hours = 1.4 ;

try

    Batt_data_LiIon_generic_singleBatt = load( 'Batt_data_LiIon_generic_singleBatt.mat' ) ;
    Batt_data_LiIon_generic_singleBatt = Batt_data_LiIon_generic_singleBatt.Batt_data ;

catch
    error( 'Remember to add complete EPOPEA repository to path' ) ;
end

Battery_data = Batt_data_LiIon_generic_singleBatt ;
Batt_capacityMargin_percent = 0 ;
powerRegulationMethod = 'DET' ;

[ TheoreticalBattCapacity_Wh_withmargin, TheoreticalBattMass_kg_withmargin, TheoreticalBattVolume_dm3_withmargin, RequiredEnergy_Wh_nomargin ] = BatterySizing_withRTG( PowerDefect_watt, PowerDefectDuration_hours, Battery_data, Batt_capacityMargin_percent, powerRegulationMethod )

SystemVoltage = 30 ;
BatteryCapacityRequired_Ah = RequiredEnergy_Wh_nomargin/SystemVoltage 
fprintf('BatteryCapacityRequired_Ah is the same as Orbilander\n')
fprintf('Mass is in the order of real batteries with this capacity in Ah\n')