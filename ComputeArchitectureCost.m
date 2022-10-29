function [f, varargout] = ComputeArchitectureCost(Archdata, costweights, varargin)
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


%% Function code

paramset = Archdata.paramset;

% DEVNOTES: 
% 1) Take vector of weights
% 2) Take vector/struct (?) of parameters values from other codes
% 3) Each parameter MUST be identifiable either by a struct field or by
% another type of ID
% 4) To include coefficient that are parameter dependent: varargin or
% struct with field "multiplier" for example


% Evaluate single cost function components
f_components = nan(length(paramset), 1);

for param_id = 1:length(paramset)
    f_components(param_id) = costweights*paramvalue; % or function to call to estimate the parameter value
end

% Evaluate total cost function
f = sum(f_components);

end