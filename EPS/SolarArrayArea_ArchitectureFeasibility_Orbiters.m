clearvars ; close all ; clc ;
set( 0, 'defaultlegendinterpreter', 'latex' ) ;
set( 0, 'defaulttextinterpreter', 'latex' ) ;
fontsize = 15 ;

color = { [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840] } ;

% LOGICAL FLOW: 
% 1) EPS_Tradeoff_SA.m script to determine which of the architectures COULD use solar arrays (evaluating an 11y EoL). 
% 2) This script considers the mission over time, only for the NSO orbiter module.

% Hypotheses:
%   1) The distance is kept constant at the Saturn distance
%   2) Assume preliminarly that there are no eclipses. This means that only
%   the FEASIBILITY of having solar arrays are evaluated. If it makes sense
%   to conduct further analysis, then this hypothesis will be removed to
%   get more refined results (once the eclipse windows are known from MA)
%   3) Solar aspect angle is 0deg for the whole lifetime, this clearly represents
%      a best case.
%   
%   NOTE: removing hypotheses 1) and 2) would only bring to larger solar arrays
%           --> If with these hypotheses use of solar arrays is
%           unrealistic, it does not make sense to further refine the
%           analysis


% Load solar array data, do not touch this, only change things in the box below
try
    SA_data_XTE_LILT_Spectrolab = load('SA_data_XTE_LILT_Spectrolab.mat') ;
    SA_data_XTE_LILT_Spectrolab = SA_data_XTE_LILT_Spectrolab.SA_data ;
    SA_data_C4MJ_CVP_Spectrolab = load('SA_data_C4MJ_CVP_Spectrolab.mat') ;
    SA_data_C4MJ_CVP_Spectrolab = SA_data_C4MJ_CVP_Spectrolab.SA_data ;
catch
    error('Remember to add complete EPOPEA repository to path')
end




% ONLY CHANGE SIZING VARIABLES IN THE BOX BELOW TO OBTAIN SOLAR ARRAY SIZE
% ------------------------------------------------------------------------

% Distances
distance_Saturn = ( 1514.50e6 + 1352.55e6 ) / 2 ; % [km] - ~9.5 AU

% Lifetime
lifetime_start = 15 ; % Estimated arrival at Saturn system                              
lifetime_end = lifetime_start+5 ; % Estimated 2.5 year moon tour + 3.5 year mission (1.5 on orbit + 2 on ground)

% Power requirement
Te = 0 ;                                          % As a first approximation, assume that during interplanetary leg you have no eclipse
Pe_watt = 0 ;                                     % Can be any number, since this is unused in case with no eclipses
Pd_watt = [ 100, 200, 400 ] ;               % [W] - Range of possible NOMINAL power requirements of each architecture, of onle the lander
Td = 1 ;                                          % Can be any number, since this is unused in case with no eclipses

% Select which solar array data to use, from loaded ones above
SA_data = SA_data_C4MJ_CVP_Spectrolab ;

% Select incidence angles to produce plots (families of curves parametrized by this)
alpha_incidence_degrees_vect = [ 0 ] ; % Assuming that in the best-case, conditions, you have always a perfect incidence of solar panel surface to Sun

% Power regulation method
powerRegulationMethod = 'DET' ; % Assume direct energy transfer for long-lifetime missions

% ------------------------------------------------------------------------
% ONLY CHANGE SIZING VARIABLES IN THE BOX ABOVE TO OBTAIN SOLAR ARRAY SIZE















% -----------------------------------------------------------------
% DO NOT TOUCH BELOW HERE :(
% -----------------------------------------------------------------

lifetime_points = 2*ceil(lifetime_end-lifetime_start)+1 ;                            % Mesh grid representing how many points since start and end of lifetime AT SATURN

lifetime = linspace( lifetime_start, lifetime_end, lifetime_points ) ; 

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

        plot( lifetime, A_SA_theoretical_SaturnSystem( : , Pd_iter, alpha_iter ), 'o-', 'linewidth', 2 ) ;
        lgtxt{Pd_iter} = [ '\textbf{NSO power requirement ', num2str(Pd_watt(Pd_iter)), ' W}' ] ;

    end
    title( ['\textbf{Sun aspect angle:} \boldmath{$\theta = ', num2str(alpha_incidence_degrees_vect(alpha_iter)), '^o$}' ], 'fontsize', fontsize ) ;
    ylabel( '\textbf{Required SA area [}\boldmath{$m^2$}\textbf{]}', 'fontsize', fontsize ) ;
    xlabel( '\textbf{Lifetime once Saturn system is reached [years]}', 'fontsize', fontsize ) ;
    ylim([0, max(A_SA_theoretical_SaturnSystem( : , Pd_iter, alpha_iter ))+20]) ;
    legend( lgtxt, 'location', 'best', 'fontsize', fontsize ) ;

end