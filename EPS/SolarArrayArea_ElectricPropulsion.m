clearvars ; close all ; clc ;
set( 0, 'defaultlegendinterpreter', 'latex' ) ;
set( 0, 'defaulttextinterpreter', 'latex' ) ;

color = { [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840] } ;

%% ELECTRIC PROPULSION: INTERPLANETARY LEG EVALUATION


% ONLY CHANGE SIZING VARIABLES IN THE BOX BELOW TO OBTAIN SOLAR ARRAY SIZE
% ------------------------------------------------------------------------

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

% ------------------------------------------------------------------------
% ONLY CHANGE SIZING VARIABLES IN THE BOX ABOVE TO OBTAIN SOLAR ARRAY SIZE















% -----------------------------------------------------------------
% DO NOT TOUCH BELOW HERE :(
% -----------------------------------------------------------------


distance_km = linspace( distance_min_km, distance_max_km, distancePoints ) ; % Create vector of distances
lifetime_years_cruise = linspace( lifetime_mindist, lifetime_maxdist, distancePoints ) ;
lifetime_years_SaturnSystem = linspace( lifetime_maxdist, lifetime_maxdist + 6, distancePoints ) ;
A_SA_theoretical_cruise = zeros( distancePoints, length(Pd_watt) ) ;
A_SA_theoretical_SaturnSystem = zeros() ;
% Outer loops spans the range given by Pd_watt, inner loop spans mission lifetimes
for j = 1:length(Pd_watt)

    for k = 1:distancePoints

        % CRUISE PHASE
        A_SA_theoretical_cruise(k,j)  = SAsizing_theoretical( distance_km(k), Pe_watt, Te, Pd_watt(j), Td, SA_data, lifetime_years_cruise(k), powerRegulationMethod ) ;

        % SATURN SYSTEM PHASE
        A_SA_theoretical_SaturnSystem(k,j) = SAsizing_theoretical( distance_km(end), Pe_watt, Te, Pd_watt(j), Td, SA_data, lifetime_years_SaturnSystem(k), powerRegulationMethod ) ;

    end

end


% Plot Area of solar array in function of distance
figure() ;
hold on ; grid on ;
lgtxt = cell(1,length(Pd_watt)) ;
for j = 1:length(Pd_watt)

    plot( lifetime_years_cruise, A_SA_theoretical_cruise(:,j), '-', 'linewidth', 1.5, 'color', color{j} ) ;
    plot( lifetime_years_SaturnSystem, A_SA_theoretical_SaturnSystem(:,j), '-', 'linewidth', 1.5, 'color', color{j}, 'handlevisibility', 'off' ) ;
    lgtxt{j} = [ '\textbf{PS power requirement: ' , num2str(Pd_watt(j)/1e3), ' [kW]}' ] ;

end
indx = floor(distancePoints/2) ;
plot( lifetime_years_cruise(end)*[ 1, 1 ], [ 0, A_SA_theoretical_SaturnSystem(end,end) ], 'k--', 'linewidth', 2 ) ;
% plot( [lifetime_years_cruise(indx)*ones(1,length(Pd_watt)), lifetime_years_SaturnSystem(1)*ones(1,length(Pd_watt))], [A_SA_theoretical_cruise(indx, 1:length(Pd_watt)), A_SA_theoretical_SaturnSystem(1,1:length(Pd_watt))], 'ko', 'markersize', 5, 'markerfacecolor', 'k' ) ;

xlabel('\textbf{Years after launch}', 'interpreter', 'latex', 'fontsize', 15 ) ;
ylabel('\textbf{Required solar array area} \boldmath{$[m^2]$}', 'interpreter', 'latex', 'fontsize', 15 ) ;
xlim([lifetime_years_cruise(1),lifetime_years_SaturnSystem(end)]) ;
text( 2, 2500, '\textbf{Cruise phase}', 'fontsize', 15 ) ;
text( 8, 2500, '\textbf{Saturn system phase}', 'fontsize', 15 ) ;
legend(lgtxt, 'fontsize', 15) ;



