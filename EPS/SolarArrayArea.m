clearvars ; close all ; clc ;

% ONLY CHANGE VARIABLES IN THE BOX BELOW TO OBTAIN SOLAR ARRAY SIZE
% -----------------------------------------------------------------


% -----------------------------------------------------------------
% ONLY CHANGE VARIABLES IN THE BOX ABOVE TO OBTAIN SOLAR ARRAY SIZE

%% ELECTRIC PROPULSION: INTERPLANETARY LEG
%       --> This section assumes you have no eclipse conditions during
%           heliocentric interplanetary leg.


%% LOCAL FUNCTIONS
% -----------------------------------------------------------------

function [] = SAsizing_interplanetary( meanDistance_km, Pe_watt, Te_seconds, Pd_watt, Td_seconds, SA_data, lifetime_years, powerRegulationMethod )

% INPUTS:
%   meanDistance_km is a number indicating the mean distance in KM of the S/C from the sun
%   Pe_watt is a number indicating the power required in eclipse
%   Te_seconds is a number indicating the duration in SECONDS of eclipse
%   Pd_watt is a number indicating the power required in daylight
%   Td_seconds is a number indicating the duration in SECONDS of daylight
%   SA_data is a structure containing data relative to the chosen solar array
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

% Unpack struct
eta_SA = SA_data.eta_SA ;
alpha_incidence_degrees = SA_data.alpha_incidence_radians ;
Id = SA_data.Id ;
yearlyDegradation_percent = SA_data.yearlyDegradation_percent ;

% Define constants
AU = 149597870.7 ; % [ km ] - Astronomical unit
P0_1AU = 1367 ; % [ W/m^2 ] - Sun irradiated power @ 1AU

% Compute Sun irradiated power @ distance of  meanDistance_km
P_SunIrradiated = P0_1AU * ( AU / meanDistance_km )^2 ;

% Compute power that the solar array can generate per m^2, at BOL
P_SA_BOL = eta_SA * P_SunIrradiated * Id * cosd( alpha_incidence_degrees ) ;

% Compute power that the solar array can generate per m^2, at EOL
P_SA_EOL = P_SA_BOL * ( 1 - yearlyDegradation_percent )^lifetime_years ; % [W/m^2]

% Compute theoretical area of solar array that is required (not considering geometrical properties of cells)
A_SA_theoretical = P_SA / P_SA_EOL ; % [m^2]

switch powerRegulationMethod

    case 'PPT'

        error('PPT has not yet been implemented, as it is better for shorter mission durations')

    case 'DET'

        Xd = 0.85 ; % Power regulation efficiency factor in daylight for DET
        Xe = 0.65 ; % Power regulation efficiency factor in eclipse for DET

end



end
