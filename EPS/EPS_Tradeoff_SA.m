clearvars ; close all ; clc ;
set( 0, 'defaultlegendinterpreter', 'latex' ) ;
set( 0, 'defaulttextinterpreter', 'latex' ) ;

color = { [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840] } ;

warning( 'Run the first section, and then run each section of interest' ) ;

% Load SA data
try
    SA_data_XTE_LILT_Spectrolab = load('SA_data_XTE_LILT_Spectrolab.mat') ;
    SA_data_XTE_LILT_Spectrolab = SA_data_XTE_LILT_Spectrolab.SA_data ;
    SA_data_C4MJ_CVP_Spectrolab = load('SA_data_C4MJ_CVP_Spectrolab.mat') ;
    SA_data_C4MJ_CVP_Spectrolab = SA_data_C4MJ_CVP_Spectrolab.SA_data ;
catch
    error('Remember to add complete EPOPEA repository to path')
end

% Define constants
SaturnMeanDistance_km = 9.6 * ASMADconstants('AU') ; % [km]

% Define system design parameters
powerRegulationMethod = 'DET' ;

%% Compute first SA sizing for NSOSL (orbiter) architecture, with no eclipses
% --------------------------------------

% Worst-case requirements
Pe_watt_NSOSL_orb = 0 ;
Te_NSOSL_orb = 0 ;
Pd_watt_NSOSL_orb = 350 ;
Td_NSOSL_orb = 1 ;
lifetime_years_NSOSL_orb = 11 ; % EoL conditions

% Size SA
SA_data = SA_data_C4MJ_CVP_Spectrolab ; % Considering future cell technologies using concentrators --> This imposes stricter pointing requirements on ADCS, with respect to no concentrators
SA_data.alpha_incidence_degrees = 0 ; % Best case
A_SA_theoretical_NSOSL_orb = SAsizing_theoretical( SaturnMeanDistance_km, Pe_watt_NSOSL_orb, Te_NSOSL_orb, Pd_watt_NSOSL_orb, Td_NSOSL_orb, SA_data, lifetime_years_NSOSL_orb, powerRegulationMethod ) ; 

%% Compute first SA sizing for NSOSL (lander) architecture, with no eclipses
% --------------------------------------

% Worst-case requirements
Pe_watt_NSOSL_lan = 0 ;
Te_NSOSL_lan = 0 ;
Pd_watt_NSOSL_lan = 800 ;
Td_NSOSL_lan = 1 ;
lifetime_years_NSOSL_lan = 11 ; % EoL conditions

SA_data = SA_data_C4MJ_CVP_Spectrolab ; % Considering future cell technologies using concentrators --> This imposes stricter pointing requirements on ADCS, with respect to no concentrators
SA_data.alpha_incidence_degrees = 0 ; % Best case
A_SA_theoretical_NSOSL_lan = SAsizing_theoretical( SaturnMeanDistance_km, Pe_watt_NSOSL_lan, Te_NSOSL_lan, Pd_watt_NSOSL_lan, Td_NSOSL_lan, SA_data, lifetime_years_NSOSL_lan, powerRegulationMethod ) ; 


%% Compute first SA sizing for SO + SL (orbiter) architecture, with no eclipses 
% ------------------------------------------

% Worst-case requirements
Pe_watt_SOSL_orb = 0 ;
Te_SOSL_orb = 0 ;
Pd_watt_SOSL_orb = 850 ;
Td_SOSL_orb = 1 ;
lifetime_years_SOSL_orb = 11 ; % EoL conditions

SA_data = SA_data_C4MJ_CVP_Spectrolab ; % Considering future cell technologies using concentrators --> This imposes stricter pointing requirements on ADCS, with respect to no concentrators
SA_data.alpha_incidence_degrees = 0 ; % Best case
A_SA_theoretical_SOSL_orb = SAsizing_theoretical( SaturnMeanDistance_km, Pe_watt_SOSL_orb, Te_SOSL_orb, Pd_watt_SOSL_orb, Td_SOSL_orb, SA_data, lifetime_years_SOSL_orb, powerRegulationMethod ) ; 

%% Compute first SA sizing for SO + SL (lander) architecture, with no eclipses 
% ------------------------------------------

% Worst-case requirements
Pe_watt_SOSL_lan = 0 ;
Te_SOSL_lan = 0 ;
Pd_watt_SOSL_lan = 800 ;
Td_SOSL_lan = 1 ;
lifetime_years_SOSL_lan = 11 ; % EoL conditions

SA_data = SA_data_C4MJ_CVP_Spectrolab ; % Considering future cell technologies using concentrators --> This imposes stricter pointing requirements on ADCS, with respect to no concentrators
SA_data.alpha_incidence_degrees = 0 ; % Best case
A_SA_theoretical_SOSL_lan = SAsizing_theoretical( SaturnMeanDistance_km, Pe_watt_SOSL_lan, Te_SOSL_lan, Pd_watt_SOSL_lan, Td_SOSL_lan, SA_data, lifetime_years_SOSL_lan, powerRegulationMethod ) ; 

