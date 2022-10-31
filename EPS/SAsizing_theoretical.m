function A_SA_theoretical = SAsizing_theoretical( meanDistance_km, Pe_watt, Te, Pd_watt, Td, SA_data, lifetime_years, powerRegulationMethod )

%% PROTOTYPE
% A_SA_theoretical = SAsizing_theoretical( meanDistance_km, Pe_watt, Te, Pd_watt, Td, SA_data, lifetime_years, powerRegulationMethod )
%
%% DESCRIPTION
% Computes required solar array area
%
%
%% INPUTS:
%   meanDistance_km: is a number indicating the mean distance in KM of the S/C from the sun
%   Pe_watt: is a number indicating the power required in eclipse
%   Te: is a number indicating the duration of the eclipse -> units not important as long as they are the same as Td
%   Pd_watt: is a number indicating the power required in daylight
%   Td: is a number indicating the duration of daylight -> units not important as long as they are the same as Te
%   SA_data: is a structure containing data relative to the chosen solar array
%       Struct contains:
%                       1)  SA_data.Id - inherent degradation factor
%                           (typically ~0.77), representing percentage of
%                           total area that actually generates power, while
%                           the remainder of the area contains
%                           interconnections, cables, etc..., between cells
%                       
%                       2) SA_data.eta_SA - Conversion efficiency of the solar
%                          cells (typically ~0.3 for triple junction GaAs, ~0.15 for monocristalline Si)
%                          
%                       3) SA_data.alpha_incidence_degrees - Mean solar incidence on solar panel
%
%                       4) SA_data.yearlyDegradation_percent - Percentage
%                          of yearly degradation of solar panels due to
%                          radiation, etc... (typically 3.75% for GaAs and 2.5% for Si)
%
%   lifetime_years is a number indicating after how many years you have the end of life conditions
%   powerRegulationMethod can be string containing 'DET' or 'PPT'

%% OUTPUT
% A_SA_theoretical: theoretical solar array area (not considering geometry of real cells)
%
%
%% CHANGELOG
% Date, User, brief summary of the modification
%
%
%% DEPENDENCIES
%
%
%% Future upgrades


% Unpack struct and define variables based on inputs
eta_SA = SA_data.eta_SA ;
alpha_incidence_degrees = SA_data.alpha_incidence_degrees ;
Id = SA_data.Id ;
yearlyDegradation_percent = SA_data.yearlyDegradation_percent ;

switch powerRegulationMethod
    case 'PPT'
        error('PPT has not yet been implemented, as it is better for shorter mission durations')
    case 'DET'
        Xd = 0.85 ; % Power regulation efficiency factor in daylight for DET
        Xe = 0.65 ; % Power regulation efficiency factor in eclipse for DET
end

% Define constants
AU = 149597870.7 ; % [ km ] - Astronomical unit
P0_1AU = 1367 ; % [ W/m^2 ] - Sun irradiated power @ 1AU

% Compute Sun irradiated power @ distance of  meanDistance_km
P_SunIrradiated = P0_1AU * ( AU / meanDistance_km )^2 ;

% Compute area required from solar array, considering daylight and eclipse conditions
P_SA = ( (Pe_watt*Te/Xe) + (Pd_watt*Td/Xd) ) / Td ; % [W]

% Compute power that the solar array can generate per m^2, at BOL
P_SA_BOL = eta_SA * P_SunIrradiated * Id * cosd( alpha_incidence_degrees ) ;

% Compute power that the solar array can generate per m^2, at EOL
P_SA_EOL = P_SA_BOL * ( 1 - yearlyDegradation_percent )^lifetime_years ; % [W/m^2]

% Compute theoretical area of solar array that is required (not considering geometrical properties of cells)
A_SA_theoretical = P_SA / P_SA_EOL ; % [m^2]

end

