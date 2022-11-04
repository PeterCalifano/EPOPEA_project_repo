function [max_turn_angle, hyp_params] = DetermineMaxTurnAngle(v_inf, bodyname, Rp_min)
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

mu = cspice_bodvrd(bodyname, 'GM', 1);
% v_inf_3rdGA = 9; % [km/s]
% Rp_mi = 400; % [km]

% Compute eccentricity of the flyby hyperbola
e_hyp = 1 + (Rp_min*v_inf^2)./mu;
% Compute maximum turning angle
turn_max = @(v_inf) 2*asind(1./e_hyp);
max_turn_angle = turn_max(v_inf);

hyp_params.e_hyp = e_hyp;

end