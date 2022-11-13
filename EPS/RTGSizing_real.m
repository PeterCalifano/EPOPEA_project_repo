function [ Number_RTGs, M_tot_RTGs, P_dissipatedThermalTotal, NuclearFuelMassTOTAL ] = RTGSizing_real( P_required_EoL, t_EoL_years, RTG_data )
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
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