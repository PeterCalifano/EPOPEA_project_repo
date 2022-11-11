clearvars ; close all ; clc ;
set( 0, 'defaultlegendinterpreter', 'latex' ) ;
set( 0, 'defaulttextinterpreter', 'latex' ) ;

color = { [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840] } ;

% Only orbiters are evaluated for now, landers will be considered
% separately as they are more critical (i.e. you need to consider visibility windows as well)
% As has been found in the electric primary propulsion sizing, the main
% criticality is present at Saturn's distance from the Sun

% Hypotheses:
%   1) Function considers that lifetime starts at 7 years, and goes on for
%   a further 6 years (about twice of what is required, but depending on
%   the Saturn Moon tour, an extra 2.5 years of lifetime could be
%   requested)
%   2) The distance is kept constant at the Saturn distance
%   3) Assume initially that there are no eclipses. This means that only
%   the FEASIBILITY of having solar arrays are evaluated. If it makes sense
%   to conduct further analysis, then this hypothesis will be removed to
%   get more refined results (once the eclipse windows are known from MA).




% ONLY CHANGE SIZING VARIABLES IN THE BOX BELOW TO OBTAIN SOLAR ARRAY SIZE
% ------------------------------------------------------------------------

warning( 'Ask Lavagna if we decide to open up Solar panels only once we arrive at Saturn, do they degrade over the previous years? ' ) ;

% Distances
distance_Saturn = ( 1514.50e6 + 1352.55e6 ) / 2 ; % [km]

% Lifetime
lifetime_start = 7 ;                              
lifetime_end = lifetime_start+6 ;                 
lifetime_points = 30 ;                            % Mesh grid representing how many points since start and end of lifetime AT SATURN

% Power requirement
Te = 0 ;                                          % As a first approximation, assume that during interplanetary leg you have no eclipse
Pe_watt = 0 ;                                     % Can be any number, since this is unused in case with no eclipses
Pd_watt = [ 300, 500, 800, 1000 ] ;               % [W] - Range of possible NOMINAL power requirements of each architecture, of onle the lander
Td = 1 ;                                          % Can be any number, since this is unused in case with no eclipses

% Solar array data
load('SA_data_XTE_LILT_Spectrolab.mat') ;
load('SA_data_C4MJ_CVP_Spectrolab.mat') ;

alpha_incidence_degrees_vect = [ 0, 7.5, 15] ; % Assuming that in the best-case, conditions, you have always a perfect incidence of solar panel surface to Sun


% Power regulation method
powerRegulationMethod = 'DET' ; % Assume direct energy transfer for long-lifetime missions

% ------------------------------------------------------------------------
% ONLY CHANGE SIZING VARIABLES IN THE BOX ABOVE TO OBTAIN SOLAR ARRAY SIZE















% -----------------------------------------------------------------
% DO NOT TOUCH BELOW HERE :(
% -----------------------------------------------------------------

lifetime = linspace( lifetime_start+2, lifetime_end, lifetime_points ) ; 

% Initialize variables for for loop
A_SA_theoretical_SaturnSystem = zeros( length(alpha_incidence_degrees_vect), length(Pd_watt), lifetime_points ) ;

for alpha_iter = 1:length(alpha_incidence_degrees_vect)

    for Pd_iter = 1:length(Pd_watt)

        for lifetime_iter = 1:length(lifetime)

            SA_data.alpha_incidence_degrees = alpha_incidence_degrees_vect(alpha_iter) ; 
            % SATURN SYSTEM PHASE
            A_SA_theoretical_SaturnSystem( lifetime_iter, Pd_iter, alpha_iter ) = SAsizing_theoretical( distance_Saturn, Pe_watt, Te, Pd_watt(Pd_iter), Td, SA_data, lifetime(lifetime_iter), powerRegulationMethod ) ;

        end

    end

end


% Plot Area of solar array in function of distance
figure() ;

for alpha_iter = 1:length(alpha_incidence_degrees_vect)

    subplot( length( alpha_incidence_degrees_vect ) , 1 , alpha_iter  ) ;
    hold on ; grid on ;
    for Pd_iter = 1:length(Pd_watt)

        plot( lifetime, A_SA_theoretical_SaturnSystem( : , Pd_iter, alpha_iter ) ) ;
        lgtxt{Pd_iter} = [ '\textbf{Power requirement ', num2str(Pd_watt(Pd_iter)), ' W}' ] ;

    end
    title( ['\textbf{Sun aspect angle:} \boldmath{$\theta = ', num2str(alpha_incidence_degrees_vect(alpha_iter)), '^o$}' ] ) ;
    ylabel( '\textbf{Required SA area [}\boldmath{$m^2$}\textbf{]}' ) ;
    xlabel( '\textbf{Lifetime once Saturn is reached}' ) ;

    legend( lgtxt, 'location', 'best' ) ;

end