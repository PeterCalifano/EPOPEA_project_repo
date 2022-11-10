function A_SA_theoretical = SAsizing_theoretical( meanDistance_km, Pe_watt, Te, Pd_watt, Td, SA_data, lifetime_years, powerRegulationMethod )

%% PROTOTYPE
% A_SA_theoretical = SAsizing_theoretical( meanDistance_km, Pe_watt, Te, Pd_watt, Td, SA_data, lifetime_years, powerRegulationMethod )
%
%% DESCRIPTION
% Computes required solar array area. Function has been validated with
% venus express mission and NASA Dawn missions, and returns quite similar
% values of solar panel area
%
%% INPUTS:
%   meanDistance_km [1,1]    Indicates mean distance in KM if the S/C with
%                            respect to the sun.
%   Pe_watt         [1,1]    Indicates the power required during eclipse
%   Te              [1,1]    Indicates the duration of the eclipse --> Units
%                            are not important as long as they are the same as Td
%   Pd_watt         [1,1]    Indicates the power required during daylight
%   Td              [1,1]    Indicates duration of daylight --> Units are not
%                            important as long as they are the same as Te
%   SA_data         [struct] Contains data of chosen solar array cells
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
%   lifetime_years        [1,1]    Indicates after how many years you have
%                                  the end of life conditions
%   powerRegulationMethod [string] String containing 'DET' or 'PPT'
%
%% OUTPUT
%   A_SA_theoretical      [1,1] theoretical solar array area (not 
%                               considering geometry of real cells)  
%
%% CHANGELOG
%   1/11/2022, Matteo D'Ambrosio, Created function
%
%% DEPENDENCIES
%   None
%
%% Future upgrades
%   Load constants from constants file
%

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

% Retrieve constants
AU = ASMADconstants('AU') ; % [km]
P0_1AU = ASMADconstants('SolarIrradiance_1AU') ; % [W/m^2]

% Compute Sun irradiated power @ distance of  meanDistance_km
P_SunIrradiated = P0_1AU * ( AU / meanDistance_km )^2 ;

% Compute area required from solar array, considering daylight and eclipse conditions
P_SA = ( (Pe_watt*Te/Xe) + (Pd_watt*Td/Xd) ) / Td ; % [W]

% Compute power that the solar array can generate per m^2, at BOL
P_SA_BOL = eta_SA * P_SunIrradiated * Id * cosd( alpha_incidence_degrees ) ; % [W/m^2]

% Compute power that the solar array can generate per m^2, at EOL
P_SA_EOL = P_SA_BOL * ( 1 - yearlyDegradation_percent )^lifetime_years ; % [W/m^2]

% Compute theoretical area of solar array that is required (not considering geometrical properties of cells)
A_SA_theoretical = P_SA / P_SA_EOL ; % [m^2]

end

