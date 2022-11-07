function [ Capacity_theoretical, Mass_theoretical, Volume_theoretical  ] = BatterySizing( Pe_watt, Te_hours, Battery_data, powerRegulationMethod )
%% PROTOTYPE 
% [Capacity_theoretical, Mass_theoretical, Volume_theoretical  ] = BatterySizing( Pe_watt, Te_hours, Battery_data, powerRegulationMethod )
%
%% DESCRIPTION
%       Function computes theoretical sizing for batteries
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

warning( 'Function has to be validated with SSEO sizing' ) ;

% Unpack struct and define variables based on inputs
N = Battery_data.NumberOfBatteries ; % Number of batteries
DOD = Battery_data.DOD ; % Depth of discharge [%]
Em = Battery_data.Em ; % Mass-specific energy [Wh/kg]
Ev = Battery_data.Ev ; % Volume-specific energy [Wh/dm^3]

switch powerRegulationMethod
    case 'PPT'
        error('PPT has not yet been implemented, as it is better for shorter mission durations')
    case 'DET'
        % Xd = 0.85 ; % Power regulation efficiency factor in daylight for DET
        Xe = 0.65 ; % Power regulation efficiency factor in eclipse for DET
end

% Retrieve constants
% -none

% Compute theoretical quantities
Capacity_theoretical = Pe_watt * Te_hours / ( DOD * N * Xe ) ;

if isfield( Battery_data, 'Em' )

    Mass_theoretical = Pe_watt * Te_hours / Em ; % [kg]

end

if isfield( Battery_data, 'Ev' )

    Volume_theoretical = Pe_watt * Te_hours / Ev ; % [dm^3]

end

end