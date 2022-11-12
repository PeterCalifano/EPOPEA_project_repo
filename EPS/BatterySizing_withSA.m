function [ TheoreticalBattCapacity_Wh_withmargin, TheoreticalBattMass_kg_withmargin, TheoreticalBattVolume_dm3_withmargin, RequiredEnergy_Wh_nomargin ] = BatterySizing_withSA( Pe_watt, Te_hours, Batt_capacityMargin_percent, powerRegulationMethod )
%% PROTOTYPE 
%
%% DESCRIPTION
%       Function computes theoretical sizing for batteries, considering
%       RTG's as power source
%
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% Date, User, brief summary of the modification
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades


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