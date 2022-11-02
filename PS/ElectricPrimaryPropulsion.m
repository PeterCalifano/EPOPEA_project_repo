clearvars ; close all ; clc ;

% ------------------------------------------------------------------------------------------------
%
% MARGINS TO ADD:
%
%   1) Safety factor of 1.5 on total impulse, lifetime, and number of 
%      cycles of the entire mission
%
%   2) Duty cycle: 90% shall be used
%
%   3) Volume of propellant tanks: computed with total propellant 
%      mass + 5% (at least)
%
% ------------------------------------------------------------------------------------------------
%
% MASSES TO TAKE INTO ACCOUNT:
%
%   1) Electric generator (power processing system) -> Low voltage needs to
%      be converted to very high voltage (~28V -> kV)
%
%   2) Propellant storage/feed system
%   
%   3) Electric thruster assembly
%
% ------------------------------------------------------------------------------------------------


%% Typical orders of magnitude to consider
%       --> Found by analyzing past missions as well as real electric 
%           propulsion technologies

P_required = [ 2, 8 ] * 1e3 ; % [ W ] - ONLY for the thrusters, not considering any other subsystems

% This order of magnitude has been considered as an input to the function
% Solar Array Area










