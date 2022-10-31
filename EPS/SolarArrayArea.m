clearvars ; close all ; clc ;

color = { [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840] } ;

%% ELECTRIC PROPULSION: INTERPLANETARY LEG EVALUATION


% ONLY CHANGE VARIABLES IN THE BOX BELOW TO OBTAIN SOLAR ARRAY SIZE
% -----------------------------------------------------------------

% Distances
distance_Venus = ( 108939000 + 107477000 ) / 2 ;  % [km]
distance_Saturn = ( 1514.50e6 + 1352.55e6 ) / 2 ; % [km]

distance_min_km = distance_Venus ;                % [km] - Minimum distance from Sun, considers possible venus flyby
distance_max_km = distance_Saturn ;               % [km] - Maximum distance from Sun,
distancePoints = 15 ;                             % Number of distances to consider, from distance_min_km to distance_max_km

% Lifetime
lifetime_mindist = 0.5 ;                           % Assumption: Year at which you reach the distance distance_min_km (Venus)  (very approximate, but is used to represent the trajectory given by MA)
lifetime_maxdist = 7 ;                             % Assumption: Year at which you reach the distance distance_max_km (Saturn) (very approximate, but is used to represent the trajectory given by MA)

% Power requirement
Te = 0 ;                                          % As a first approximation, assume that during interplanetary leg you have no eclipse
Pe_watt = 0 ;                                     % Can be any number, since this is unused in case with no eclipses
Pd_watt = [ 2, 4, 6, 8 ] * 1e3 ;                  % [W] - Range of possible NOMINAL power requirements, considering ONLY THE POWER REQUIRED BY THE ELECTRIC PROPULSION SYSTEM 
Td = 1 ;                                          % Can be any number, since this is unused in case with no eclipses

% Solar array data
%       --> Assuming typical values for triple-junction GaAs solar panels
%       --> Check function SAsizing_theoretical.m for description of each variable
SA_data.Id = 0.77 ; % Inherent degradation
SA_data.eta_SA = 0.3 ; % Solar panel conversion efficiency
SA_data.alpha_incidence_degrees = 0 ; % Assuming that in the best-case, conditions, you have always a perfect incidence of solar panel surface to Sun
SA_data.yearlyDegradation_percent = 0.0375 ; % Yearly degradation percentage of solar panel efficiency


% Power regulation method
powerRegulationMethod = 'DET' ; % Assume direct energy transfer for long-lifetime missions

% -----------------------------------------------------------------
% ONLY CHANGE VARIABLES IN THE BOX ABOVE TO OBTAIN SOLAR ARRAY SIZE




% -----------------------------------------------------------------
% DO NOT TOUCH BELOW HERE :(
% -----------------------------------------------------------------


distance_km = linspace( distance_min_km, distance_max_km, distancePoints ) ; % Create vector of distances
lifetime_years = linspace( lifetime_mindist, lifetime_maxdist, distancePoints ) ;
A_SA_theoretical = zeros( distancePoints, length(Pd_watt) ) ;

% Outer loops spans the range given by Pd_watt, inner loop spans mission lifetimes
for j = 1:length(Pd_watt)

    for k = 1:distancePoints

        A_SA_theoretical(k,j) = SAsizing_theoretical( distance_km(k), Pe_watt, Te, Pd_watt(j), Td, SA_data, lifetime_years(k), powerRegulationMethod ) ;

    end

end

% Plot Area of solar array in function of distance
figure() ;
hold on ; grid on ;
lgtxt = cell(1,length(Pd_watt)) ;
for j = 1:length(Pd_watt)

    plot( lifetime_years, A_SA_theoretical(:,j), '-', 'linewidth', 1.5 ) ;
    lgtxt{j} = [ 'PS power requirement: ' , num2str(Pd_watt(j)/1e3), ' [kW]' ] ;

end
xlabel('\textbf{Years after launch}', 'interpreter', 'latex', 'fontsize', 15 ) ;
ylabel('\textbf{Solar array area} \boldmath{$[m^2]$}', 'interpreter', 'latex', 'fontsize', 15 ) ;
legend(lgtxt, 'fontsize', 15) ;



%% LOCAL FUNCTIONS
% -----------------------------------------------------------------

function A_SA_theoretical = SAsizing_theoretical( meanDistance_km, Pe_watt, Te, Pd_watt, Td, SA_data, lifetime_years, powerRegulationMethod )

% INPUTS:
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




