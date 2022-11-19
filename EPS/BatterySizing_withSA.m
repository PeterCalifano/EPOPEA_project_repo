function [ TheoreticalBattCapacity_Wh_withmargin, TheoreticalBattMass_kg_withmargin, TheoreticalBattVolume_dm3_withmargin, RequiredEnergy_Wh_nomargin ] = BatterySizing_withSA( Pe_watt, Te_hours, Batt_capacityMargin_percent, powerRegulationMethod )
%% PROTOTYPE 
% [ TheoreticalBattCapacity_Wh_withmargin, TheoreticalBattMass_kg_withmargin, TheoreticalBattVolume_dm3_withmargin, RequiredEnergy_Wh_nomargin ] = BatterySizing( PowerDefect_watt, PowerDefectDuration_hours, Battery_data, Batt_capacityMargin_percent, powerRegulationMethod )
%% DESCRIPTION
%       Function computes theoretical sizing for batteries, considering
%       RTG as power source
%       NOTE: Function has been validated with Orbilander concept mission
%
%% INPUT
% Pe_watt [1x1] Worst-case power that is required during eclipse times
% Te_hours [1x1] Duration of eclipse, in worst-case conditions
% Batt_capacityMargin_percent [1x1] Percentage that represents margin considered on battery capacity
% powerRegulationMethod [string] String containing 'DET' or 'PPT'
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% TheoreticalBattCapacity_Wh_withmargin [1x1] Sized capacity of battery in Wh, considering DOD, number of batteries, and efficiency
% TheoreticalBattMass_kg_withmargin [1x1] Mass of sized battery
% TheoreticalBattVolume_dm3_withmargin [1x1] Volume of sized battery
% RequiredEnergy_Wh_nomargin [1x1] Energy that is actually required by loads
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 12/11/22 Matteo D'Ambrosio, Created function and validated it 
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -None
% -------------------------------------------------------------------------------------------------------------

% Unpack struct and define variables based on inputs
N = Battery_data.NumberOfBatteries ; % Number of batteries
DOD = Battery_data.DOD ; % Depth of discharge [%]
Em = Battery_data.Em ; % Mass-specific energy [Wh/kg]
Ev = Battery_data.Ev ; % Volume-specific energy [Wh/dm^3]

switch powerRegulationMethod

    case 'PPT'

        error('PPT has not yet been implemented, as it is better for shorter mission durations')

    case 'DET'
        % Xd = 0.85 ; % Power regulation efficiency factor while power source is providing power for DET
        Xe = 0.65 ; % Power regulation efficiency factor while battery is discharging for DET -> This is the one to consider while battery is powering the loads alongside the RTG

end

% Retrieve constants
% -none

% Compute required energy --> Energy that actually needs to be provided to
% the loads
RequiredEnergy_Wh_nomargin = Pe_watt * Te_hours ;

% Compute battery capacity with and without margin
TheoreticalBattCapacity_Wh_nomargin = RequiredEnergy_Wh_nomargin / ( DOD * N * Xe ) ; % [Wh]
TheoreticalBattCapacity_Wh_withmargin = TheoreticalBattCapacity_Wh_nomargin * ( 1 + Batt_capacityMargin_percent ) ; % [Wh]

% Compute theoretical mass and volume (not considering size and shape of battery cells)
if isfield( Battery_data, 'Em' )

    TheoreticalBattMass_kg_withmargin = TheoreticalBattCapacity_Wh_withmargin / Em ; % [kg]

end

if isfield( Battery_data, 'Ev' )

    TheoreticalBattVolume_dm3_withmargin = TheoreticalBattCapacity_Wh_withmargin / Ev ; % [dm^3]

end

end