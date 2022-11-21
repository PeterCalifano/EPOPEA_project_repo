function [nSeries, nParallel, V_real, C_real] = refinedBatterySizing_withRTG(Cell_data, C_required, V_required)
%% PROTOTYPE 
% [nSeries, nParallel, V_real, C_real] = refinedBatterySizing_withRTG(Cell_data, C_required, V_required)
%% DESCRIPTION
%       Function computes refined sizing for batteries, considering
%       RTG as power source and data coming from the output of the
%       theoretical sizing of batteries, see function
%       BatterySizing_withRTG.
%
%% INPUT
% Cell_data [struct]: contains the useful data for a single cell
% C_required [1x1]: fequired capacity of the battery, as for output of the
% BatterySizing_withRTG function (it's the total energy in Wh)
% V_required [1x1]: required system voltage. Depends on the SC needs.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% nSeries [1x1]: number of cells in a string
% nParallel [1x1]: number of strings in parallel
% V_real [1x1]: real voltage of the battery (as voltage of a single string)
% C_real [1x1]: real capacity of the battery (sum of the string capacities, in Wh)
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19/11/22 Davide Intra, Created function  
% -------------------------------------------------------------------------------------------------------------

% Data unpacking
V_cell = Cell_data.V; % Voltage [V]
C_cell = Cell_data.C; % Capacity [Ah]

% Define number of cells per string
nSeries = ceil(V_required/V_cell);

% Define real voltage of the string
V_real = nSeries*V_cell;

% Define capacity of ONE string
mu = 0.8; % packing efficiency
C_string = mu*C_cell*V_real; % [Wh]

% Define number of series
nParallel = ceil(C_required/C_string);

% Define real battery capacity
C_real = nParallel*C_string; % [Wh]

% Check in case of one cell failure -> one string is lost
if C_real/nParallel * (nParallel-1) - C_required < 0
    nParallel = nParallel+1; % add one string to be sure
end

% Update new C_real
C_real = nParallel*C_string;

end

