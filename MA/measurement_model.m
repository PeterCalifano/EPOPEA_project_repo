function [range, azimuth, elevation] = measurement_model(station_name,time,rv_eci)
% 
%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The function implements the measurement model to obtain the measures
%   in the topocentric frame of a given station from the position of the
%   object in inertial frame at a certain time.
%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% station_name - Name of the station 
% time - Epoch at which the measurement happens [s]
% rv_eci - State of the object in ECI [Km, Km/s]
%
%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% range - Range with respect to the ground station topocentric frame [Km]
% azimuth - Azimuth with respect to the ground station topocentric frame [deg]
% elevation - Elevation with respect to the ground station topocentric frame [deg]
% 
%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matteo Lusvarghi
%
% -------------------------------------------------------------------------

% Compute the position of the station in ECI
rv_station_eci = cspice_spkezr(station_name, time, 'J2000', 'NONE', 'EARTH');

% Compute relative position of the object wrt the station in ECI
rv_rel_eci = rv_eci - rv_station_eci;

% Compute rotation matrix from ECI to topo frame of the station
ROT_ECI2TOPO = cspice_sxform('J2000', [station_name,'_TOPO'] , time);

% Compute relative position of the object wrt the station in TOPO
rv_rel_topo = ROT_ECI2TOPO * rv_rel_eci;

% Compute range, Azimuth and Elevation
[range, azimuth,  elevation] = cspice_reclat(rv_rel_topo(1:3));

end