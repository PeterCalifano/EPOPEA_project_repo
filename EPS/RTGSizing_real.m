function [ Number_RTGs, M_tot_RTGs, P_dissipatedThermalTotal, NuclearFuelMassTOTAL ] = RTGSizing_real( P_required_EoL, t_EoL_years, RTG_data )
%% PROTOTYPE
% [ Number_RTGs, M_tot_RTGs, P_dissipatedThermalTotal, NuclearFuelMassTOTAL ] = RTGSizing_real( P_required_EoL, t_EoL_years, RTG_data )
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computes number, mass and dissipated power of specific RTG necessary to obtain Power P_required_EoL at time t_EoL_years.
% NOTE: Function has been validated with Mars Perseverance RTG's 
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% P_required_EoL [1x1] Amount of power that is required by the loads at EoL
% t_EoL_years [1x1] Years after start of mission, at which power P_required_EoL is requested (e.g. EoL)
% RTG_data [struct] Structure containing information on RTG data. This structure is obtained by loading a variable from RTG_data folder
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% Number_RTGs [1x1] Number of RTGs that are necessary
% M_tot_RTGs [1x1] Total mass of RTGs that are required on mission
% P_dissipatedThermalTotal [1x1] Total thermal power that is dissipated by the RTG during energy conversion
% NuclearFuelMassTOTAL [1x1] Total mass of nuclear fuel that is contained in RTGs
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 13/11/2022 Matteo D'Ambrosio, Created and validated function
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -None
% -------------------------------------------------------------------------------------------------------------


% Unpack structure - "unit" stands for a single RTG unit
P_BoL_electric_unit = RTG_data.BOLpower_electric ;
P_BoL_thermal_unit = RTG_data.BOLpower_thermal ;
M_RTG_unit = RTG_data.RTGmass ;
Fuel_mass_unit = RTG_data.FUELmass ;
t_halflife_years = RTG_data.halflife ;

% Compute EoL power output of unit
powerReductionFactor = 2^(-t_EoL_years/t_halflife_years) ; % Formula that represents exponential power decay of isotope
P_EoL_electric_unit = P_BoL_electric_unit * powerReductionFactor ;
P_EoL_thermal_unit = P_BoL_thermal_unit * powerReductionFactor ;

% Compute number of RTGs required
Number_RTGs = ceil( P_required_EoL / P_EoL_electric_unit ) ;

% Compute total mass of units required
M_tot_RTGs = Number_RTGs * M_RTG_unit ;

% Compute EoL dissipated thermal power
P_dissipatedThermalTotal = P_EoL_thermal_unit * Number_RTGs ;

% Compute total nuclear fuel mass required
NuclearFuelMassTOTAL = Number_RTGs * Fuel_mass_unit ;

end