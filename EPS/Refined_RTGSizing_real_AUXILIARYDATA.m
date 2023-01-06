function [ P_EoL_electric_total, M_tot_RTGs, P_dissipatedThermalTotal, NuclearFuelMassTOTAL ] = Refined_RTGSizing_real_AUXILIARYDATA( t_EoL_years, Number_RTGs1, Number_RTGs2, RTGdata1, RTGdata2 )
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
M_RTG_unit1 = RTGdata1.RTGmass ;
Fuel_mass_unit1 = RTGdata1.FUELmass ;
M_RTG_unit2 = RTGdata2.RTGmass ;
Fuel_mass_unit2 = RTGdata2.FUELmass ;
t_halflife_years = RTGdata1.halflife ;


% Compute BoL power output of unit
P_BoL_electric_unit1 = RTGdata1.BOLpower_electric ;
P_BoL_electric_unit2 = RTGdata2.BOLpower_electric ;
P_BoL_electric_thermal1 = RTGdata1.BOLpower_thermal ;
P_BoL_electric_thermal2 = RTGdata2.BOLpower_thermal ;

% Compute EoL power output of unit
powerReductionFactor = 2^(-t_EoL_years/t_halflife_years) ; % Formula that represents exponential power decay of isotope
P_EoL_electric_unit1 = P_BoL_electric_unit1 * powerReductionFactor ;
P_EoL_electric_unit2 = P_BoL_electric_unit2 * powerReductionFactor ;

% Compute number of RTGs required
P_EoL_electric_total = P_EoL_electric_unit1 * Number_RTGs1 + P_EoL_electric_unit2 * Number_RTGs2 ;

% Compute total mass of units required
M_tot_RTGs = Number_RTGs1 * M_RTG_unit1 + Number_RTGs2 * M_RTG_unit2 ;

% Compute EoL dissipated thermal power
P_dissipatedThermalTotal = Number_RTGs1 * P_BoL_electric_thermal1 + Number_RTGs2 * P_BoL_electric_thermal2 ;

% Compute total nuclear fuel mass required
NuclearFuelMassTOTAL = Number_RTGs1 * Fuel_mass_unit1 + Number_RTGs2 * Fuel_mass_unit2 ;

end